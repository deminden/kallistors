//! Index loading and metadata assembly for each mapping mode.

use super::*;

pub(super) fn build_bifrost_index_with_positions_threaded_impl(
    path: &Path,
    load_positional_info: bool,
    _threads: usize,
    store_ec_blocks: bool,
    packed_graph: bool,
    defer_names: bool,
) -> Result<BifrostIndex> {
    let _header_timing = timing::scoped(Stage::IndexHeaderParse);
    let file = File::open(path).map_err(|_| Error::MissingFile(path.to_path_buf()))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    let index_version = read_u64_le(&mut reader)?;
    let dbg_size = read_u64_le(&mut reader)?;
    let dbg_start = reader.stream_position()?;
    drop(_header_timing);

    let _graph_timing = timing::scoped(Stage::GraphDecode);
    let file_format_version = read_u64_le(&mut reader)?;
    let k = read_i32_le(&mut reader)?;
    let g = read_i32_le(&mut reader)?;
    let _ = file_format_version;

    if k <= 0 || g <= 0 {
        return Err(Error::InvalidFormat("invalid k/g".into()));
    }
    let k = k as usize;
    let g = g as usize;
    if k > 32 || g > 32 {
        return Err(Error::UnsupportedFeature(
            "k/g > 32 not supported in Bifrost path".into(),
        ));
    }

    let v_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut unitigs = Vec::with_capacity(if packed_graph { 0 } else { v_unitigs_sz });
    let mut compressed_scratch = Vec::new();
    let mut packed_values = Vec::new();
    let mut packed_offsets = Vec::new();
    let mut packed_lengths = Vec::new();
    for _ in 0..v_unitigs_sz {
        let len = read_u64_le(&mut reader)? as usize;
        let data_sz = len.div_ceil(4);
        if packed_graph {
            packed_offsets.push(packed_values.len());
            packed_lengths.push(len);
            append_buffered_compressed(
                &mut reader,
                len,
                &mut packed_values,
                &mut compressed_scratch,
            )?;
        } else {
            compressed_scratch.resize(data_sz, 0);
            reader.read_exact(&mut compressed_scratch)?;
            unitigs.push(decode_compressed_sequence(&compressed_scratch, len));
        }
    }
    drop(compressed_scratch);
    let packed_bases = packed_graph.then(|| {
        packed_values.shrink_to_fit();
        EncodedBaseStore {
            metadata: SequenceMetadata::new(packed_offsets, packed_lengths),
            values: packed_values,
        }
    });
    let km_unitigs_sz = read_u64_le(&mut reader)? as usize;
    let mut km_unitigs = Vec::with_capacity(if packed_graph { 0 } else { km_unitigs_sz });
    let mut packed_km_words = Vec::with_capacity(if packed_graph { km_unitigs_sz } else { 0 });
    for _ in 0..km_unitigs_sz {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        if packed_graph {
            packed_km_words.push(u64::from_le_bytes(buf));
        } else {
            km_unitigs.push(decode_kmer(&buf, k));
        }
    }

    let h_kmers_ccov_sz = read_u64_le(&mut reader)? as usize;
    let mut h_kmers = Vec::with_capacity(if packed_graph { 0 } else { h_kmers_ccov_sz });
    let mut packed_h_words = Vec::with_capacity(if packed_graph { h_kmers_ccov_sz } else { 0 });
    let mut h_kmer_map = (h_kmers_ccov_sz > 0).then(|| HashMap::with_capacity(h_kmers_ccov_sz));
    for id in 0..h_kmers_ccov_sz {
        let word = read_u64_le(&mut reader)?;
        let fwd = word.checked_shr((64 - 2 * k) as u32).unwrap_or(0);
        let rev = crate::dna::reverse_complement_word(fwd, k);
        h_kmer_map
            .as_mut()
            .expect("nonempty abundant k-mers")
            .insert(fwd.min(rev), v_unitigs_sz + km_unitigs_sz + id);
        if packed_graph {
            packed_h_words.push(word);
        } else {
            h_kmers.push(decode_kmer(&word.to_le_bytes(), k));
        }
    }

    let _bfg_header = read_u64_le(&mut reader)?;
    let _checksum = read_u64_le(&mut reader)?;
    let _v_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _km_unitigs_sz2 = read_u64_le(&mut reader)?;
    let _h_kmers_sz2 = read_u64_le(&mut reader)?;
    let hmap_min_unitigs_sz = read_u64_le(&mut reader)? as usize;

    let nb_bmp_unitigs = read_u64_le(&mut reader)? as usize;
    let mut bmp_unitigs = Vec::with_capacity(nb_bmp_unitigs);
    for _ in 0..nb_bmp_unitigs {
        bmp_unitigs.push(BitmapPositions::new(
            read_bitcontainer_values(&mut reader)?,
            packed_graph,
        ));
    }

    let nb_bmp_km = read_u64_le(&mut reader)? as usize;
    let mut bmp_km = Vec::with_capacity(nb_bmp_km);
    for _ in 0..nb_bmp_km {
        bmp_km.push(BitmapPositions::new(
            read_bitcontainer_values(&mut reader)?,
            false,
        ));
    }

    let nb_special_minz = read_u64_le(&mut reader)? as usize;
    let nb_bmp_special = (nb_special_minz >> 32) + 1;
    let mut abundant = HashSet::new();
    let mut overcrowded = HashSet::new();
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            abundant.insert(id as usize);
        }
    }
    for _ in 0..nb_bmp_special {
        for id in read_bitcontainer_values(&mut reader)? {
            overcrowded.insert(id as usize);
        }
    }

    let mut special_minz = Vec::with_capacity(nb_special_minz);
    for i in 0..nb_special_minz {
        let minz = read_u64_le(&mut reader)?;
        let is_abundant = abundant.contains(&i);
        let is_overcrowded = overcrowded.contains(&i);
        let pos_id = if is_abundant || is_overcrowded {
            let mut id = 0xffff_ffff_0000_0000 | if is_overcrowded { 0x8000_0000 } else { 0 };
            if is_abundant {
                let count = read_u32_le(&mut reader)? as u64;
                id |= count;
            }
            id
        } else {
            read_u64_le(&mut reader)?
        };
        special_minz.push((minz, pos_id));
    }

    drop(abundant);
    drop(overcrowded);
    reader.seek(SeekFrom::Start(dbg_start + dbg_size))?;

    let mphf_size = read_u64_le(&mut reader)? as usize;
    let mut mphf_buf = vec![0u8; mphf_size];
    reader.read_exact(&mut mphf_buf)?;
    let mut mphf_cur = std::io::Cursor::new(mphf_buf.as_slice());
    let mut mphf = BooPhf::load(&mut mphf_cur)?;
    if packed_graph || (!store_ec_blocks && !load_positional_info) {
        mphf.release_redundant_rank_storage();
    }
    crate::index::scratch_pages::retire(mphf_buf);

    let dlist_size = read_u64_le(&mut reader)?;
    let _dlist_overhang = read_u64_le(&mut reader)?;
    let mut dlist = if dlist_size > 0 {
        HashSet::with_capacity(dlist_size as usize)
    } else {
        HashSet::new()
    };
    for _ in 0..dlist_size {
        let mut buf = [0u8; 8];
        reader.read_exact(&mut buf)?;
        let seq = decode_kmer(&buf, k);
        if let Some((fwd, rev)) = encode_kmer_pair(&seq) {
            dlist.insert(fwd.min(rev));
        }
    }
    let dlist = if dlist.is_empty() { None } else { Some(dlist) };

    let node_count = read_u64_le(&mut reader)? as usize;
    let mut ec_blocks = if store_ec_blocks {
        Vec::with_capacity(node_count)
    } else {
        Vec::new()
    };
    let mut flat_ec = FlatEcIndex::new(if store_ec_blocks { 0 } else { node_count });
    let mut node_buf = Vec::new();
    let abundant_base = v_unitigs_sz + km_unitigs_sz;
    if node_count != abundant_base + h_kmers_ccov_sz {
        return Err(Error::InvalidFormat(
            "unitig and node counts mismatch".into(),
        ));
    }
    let mut node_kmer = vec![0; k];
    for node_row in 0..node_count {
        if node_row < abundant_base {
            // Vector-backed graph entries retain their serialized order.
            reader.seek_relative(k as i64)?;
        } else {
            // Bifrost reinserts abundant k-mers into a hash table when loading
            // the graph. Its iteration order can differ from the graph file's
            // order. Node records carry the authoritative k-mer/EC association.
            reader.read_exact(&mut node_kmer)?;
            if packed_graph {
                associate_abundant_word(
                    &mut packed_h_words,
                    h_kmer_map.as_mut().expect("nonempty abundant k-mers"),
                    abundant_base,
                    node_row,
                    &node_kmer,
                )?;
            } else {
                associate_abundant_node(
                    &mut h_kmers,
                    h_kmer_map.as_mut().expect("nonempty abundant k-mers"),
                    abundant_base,
                    node_row,
                    &node_kmer,
                )?;
            }
        }
        let node_size = read_u32_le(&mut reader)? as usize;
        if !store_ec_blocks && !load_positional_info {
            flat_ec.push_buffered_node(&mut reader, node_size, &mut node_buf)?;
            continue;
        }
        node_buf.resize(node_size, 0);
        reader.read_exact(&mut node_buf)?;
        let mut cur = std::io::Cursor::new(node_buf.as_slice());
        let _node_id = read_u32_le(&mut cur)?;
        let blocks = if load_positional_info {
            crate::index::read_block_array_blocks_with_positions(&mut cur)?
        } else {
            read_block_array_blocks(&mut cur)?
        };
        if store_ec_blocks {
            ec_blocks.push(blocks.into_boxed_slice());
        } else {
            flat_ec.push_unitig_blocks(&blocks);
        }
    }
    flat_ec.owned_ec_blocks = ec_blocks;
    crate::index::scratch_pages::retire(node_buf);
    let load_threads = _threads.max(1);
    let rayon_pool = if load_threads > 1 {
        Some(
            ThreadPoolBuilder::new()
                .num_threads(load_threads)
                .build()
                .map_err(|e| Error::InvalidFormat(format!("failed to build rayon pool: {e}")))?,
        )
    } else {
        None
    };
    let encoded_unitigs = if let Some(unitig_bases) = packed_bases {
        EncodedUnitigs {
            unitig_bases,
            km_bases: EncodedBaseStore::from_fixed_words(packed_km_words, k),
            g,
        }
    } else {
        precompute_encoded_unitigs(&unitigs, &km_unitigs, g, rayon_pool.as_ref())
    };

    let selective_fallback_kmers = if packed_graph {
        SelectiveFallbackKmerIndex::build_packed(
            v_unitigs_sz,
            &encoded_unitigs.km_bases,
            &packed_h_words,
            k,
        )
    } else {
        SelectiveFallbackKmerIndex::build(v_unitigs_sz, &km_unitigs, &h_kmers)
    };
    crate::index::scratch_pages::retire(packed_h_words);
    let use_one_pass = packed_graph
        && rayon_pool.as_ref().is_some_and(|pool| {
            one_pass::eligible(
                &one_pass::Context {
                    encoded: &encoded_unitigs,
                    mphf: &mphf,
                    pool,
                    buckets: hmap_min_unitigs_sz,
                    k,
                },
                &bmp_unitigs,
                &bmp_km,
                special_minz.len(),
            )
        });
    let (minz_positions, minz_fingerprints) = if use_one_pass {
        drop(_graph_timing);
        one_pass::build(
            one_pass::Context {
                encoded: &encoded_unitigs,
                mphf: &mphf,
                pool: rayon_pool.as_ref().expect("eligible pool"),
                buckets: hmap_min_unitigs_sz,
                k,
            },
            bmp_unitigs,
            bmp_km,
            special_minz,
        )
    } else {
        let mut minz_counts = Vec::with_capacity(hmap_min_unitigs_sz + 1);
        minz_counts.resize(hmap_min_unitigs_sz, 0);
        drop(_graph_timing);
        let _count_timing = timing::scoped(Stage::MinimizerCountPass);
        let mut unitig_bucket_cache = (bmp_unitigs.len() == 1
            && hmap_min_unitigs_sz < u32::MAX as usize
            && (packed_graph || (!store_ec_blocks && !load_positional_info)))
            .then(|| Vec::with_capacity(bmp_unitigs[0].len()));
        // The same conservative query filter benefits the serial scanner.
        // Its construction already supports the calling thread without a pool.
        let fingerprinting =
            k == 31 && g == 23 && unitig_bucket_cache.is_some() && bmp_km.len() <= 1;
        let minz_fingerprints: Vec<_> = if fingerprinting {
            (0..hmap_min_unitigs_sz)
                .map(|_| AtomicU16::new(0))
                .collect()
        } else {
            Vec::new()
        };
        populate_unitig_minimizer_counts(
            &mut minz_counts,
            &bmp_unitigs,
            &encoded_unitigs,
            &mphf,
            rayon_pool.as_ref(),
            unitig_bucket_cache.as_mut(),
            &minz_fingerprints,
        );

        let km_glen = k.saturating_sub(g) + 1;
        let kmer_load_ctx = KmerLoadContext {
            km_bases: &encoded_unitigs.km_bases,
            g,
            km_glen,
            mphf: &mphf,
        };
        populate_kmer_minimizer_counts(
            &mut minz_counts,
            &bmp_km,
            &kmer_load_ctx,
            rayon_pool.as_ref(),
            &minz_fingerprints,
        );

        let mut inline_escapes = Vec::new();
        for (minz, pos_id) in special_minz.iter().copied() {
            let bytes = minz.to_le_bytes();
            if let Some(idx) = mphf.lookup(&bytes) {
                minz_counts[idx as usize] += 1;
                if let Some(slot) = minz_fingerprints.get(idx as usize) {
                    slot.store(super::super::fingerprints::DISABLED, Ordering::Relaxed);
                }
                if pos_id >> 63 != 0 {
                    inline_escapes.push(idx as usize);
                }
            }
        }
        inline_escapes.sort_unstable();
        inline_escapes.dedup();
        let inline_enabled = (packed_graph || (!store_ec_blocks && !load_positional_info))
            && bmp_unitigs.len() <= 1
            && bmp_km.len() <= 1
            && encoded_unitigs.unitig_bases.metadata.len() < (1usize << 31)
            && encoded_unitigs.km_bases.metadata.len() < (1usize << 31)
            && hmap_min_unitigs_sz < u32::MAX as usize
            && (bmp_unitigs.is_empty() || unitig_bucket_cache.is_some() || rayon_pool.is_none());
        let mut minz_positions =
            MinzPositionIndex::from_counts_inline(minz_counts, &inline_escapes, inline_enabled);
        let mut minz_cursors = minz_positions.cursors();
        drop(_count_timing);
        let _fill_timing = timing::scoped(Stage::MinimizerFillPass);
        populate_unitig_minimizer_storage(
            &mut minz_positions,
            &mut minz_cursors,
            &bmp_unitigs,
            &encoded_unitigs,
            &mphf,
            rayon_pool.as_ref(),
            unitig_bucket_cache.as_deref(),
        );
        if let Some(cache) = unitig_bucket_cache {
            crate::index::scratch_pages::retire(cache);
        }
        for bitmap in bmp_unitigs {
            bitmap.retire();
        }
        populate_kmer_minimizer_storage(
            &mut minz_positions,
            &mut minz_cursors,
            &bmp_km,
            &kmer_load_ctx,
            rayon_pool.as_ref(),
        );

        for bitmap in bmp_km {
            bitmap.retire();
        }
        for (minz, pos_id) in special_minz {
            let bytes = minz.to_le_bytes();
            if let Some(idx) = mphf.lookup(&bytes) {
                minz_positions.push_at(&mut minz_cursors, idx as usize, pos_id);
            }
        }

        crate::index::scratch_pages::retire(minz_cursors);
        minz_positions.finish_inline();
        (minz_positions, minz_fingerprints)
    };
    let (onlist, deferred_names) = if defer_names {
        let (data, names) = super::super::onlist::read_onlist_scanned(&mut reader)?;
        (
            data,
            names
                .map(|(start, count)| {
                    reader.get_ref().try_clone().map(|file| {
                        super::super::onlist::TranscriptNames::deferred(file, start, count)
                    })
                })
                .transpose()?,
        )
    } else {
        (read_onlist_with_lengths(&mut reader)?, None)
    };

    let cache_id = loop {
        let id = NEXT_CACHE_ID.load(Ordering::Relaxed);
        let next = id
            .checked_add(1)
            .expect("index cache identity space exhausted");
        if NEXT_CACHE_ID
            .compare_exchange_weak(id, next, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            break id;
        }
    };
    Ok(BifrostIndex {
        cache_id,
        index_version,
        k,
        g,
        unitigs,
        km_unitigs,
        h_kmers,
        h_kmer_map,
        dlist,
        encoded_unitigs,
        kmer_pos_index: None,
        kmer_index: None,
        flat_ec,
        minz_positions,
        large_candidate_cache: Vec::new(),
        candidate_cache_slots: Vec::new(),
        minz_fingerprints: minz_fingerprints
            .into_iter()
            .map(AtomicU16::into_inner)
            .collect(),
        mphf,
        selective_fallback_kmers,
        onlist: onlist.onlist,
        transcript_lengths: onlist.lengths,
        transcript_names: onlist.names,
        deferred_names,
        shade_to_color: onlist.shade_to_color,
        shade_sequences: onlist.shade_sequences,
        use_shade: onlist.use_shade,
    })
}
