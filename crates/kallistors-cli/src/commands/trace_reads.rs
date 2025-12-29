use anyhow::{Result, anyhow};
use kallistors::io::ReadSource;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn normalize_header(input: &str) -> String {
    input
        .trim()
        .strip_prefix('@')
        .unwrap_or(input.trim())
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

fn format_ec(ec: &Option<Vec<u32>>) -> String {
    match ec {
        Some(values) if !values.is_empty() => values
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
        Some(_) => "-".to_string(),
        None => "-".to_string(),
    }
}

fn load_gene_map(path: &Path) -> Result<HashMap<u32, String>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(HashMap::new());
    }
    let header = header.trim_end();
    let cols: Vec<&str> = header.split('\t').collect();
    let id_idx = cols
        .iter()
        .position(|&c| c == "transcript_id")
        .ok_or_else(|| anyhow!("gene map missing 'transcript_id' column"))?;
    let name_idx = cols
        .iter()
        .position(|&c| c == "gene_name")
        .ok_or_else(|| anyhow!("gene map missing 'gene_name' column"))?;

    let mut map = HashMap::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end();
        if row.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = row.split('\t').collect();
        if let (Some(id), Some(name)) = (parts.get(id_idx), parts.get(name_idx)) {
            let tid = id.trim();
            if let Ok(parsed) = tid.parse::<u32>() {
                map.insert(parsed, name.to_string());
                line.clear();
                continue;
            }
            if let Some((num, _)) = tid.strip_prefix("ENST").and_then(|v| v.split_once('.'))
                && let Ok(parsed) = num.parse::<u32>()
            {
                map.insert(parsed, name.to_string());
            }
        }
        line.clear();
    }
    Ok(map)
}

fn format_gene_names(ec: &Option<Vec<u32>>, map: &HashMap<u32, String>) -> String {
    let Some(values) = ec else {
        return "-".to_string();
    };
    if values.is_empty() {
        return "-".to_string();
    }
    let mut names: Vec<String> = values
        .iter()
        .filter_map(|id| map.get(id).cloned())
        .collect();
    if names.is_empty() {
        return "-".to_string();
    }
    names.sort();
    names.dedup();
    names.join(",")
}

fn load_minimizer_positions(path: &Path) -> Result<HashMap<String, Vec<usize>>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(HashMap::new());
    }
    let header = header.trim_end();
    let cols: Vec<&str> = header.split('\t').collect();
    let read_idx = cols
        .iter()
        .position(|&c| c == "read")
        .ok_or_else(|| anyhow!("minimizer positions missing 'read' column"))?;
    let pos_idx = cols
        .iter()
        .position(|&c| c == "kallisto_hit_no_kallistors_hit_positions")
        .ok_or_else(|| {
            anyhow!("minimizer positions missing 'kallisto_hit_no_kallistors_hit_positions' column")
        })?;

    let mut map: HashMap<String, Vec<usize>> = HashMap::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let row = line.trim_end();
        if row.is_empty() {
            line.clear();
            continue;
        }
        let parts: Vec<&str> = row.split('\t').collect();
        let Some(read) = parts.get(read_idx) else {
            line.clear();
            continue;
        };
        let Some(pos_text) = parts.get(pos_idx) else {
            line.clear();
            continue;
        };
        if *pos_text == "-" || pos_text.is_empty() {
            line.clear();
            continue;
        }
        let key = normalize_header(read);
        let positions = pos_text
            .split(',')
            .filter_map(|v| v.parse::<usize>().ok())
            .collect::<Vec<_>>();
        if !positions.is_empty() {
            let entry = map.entry(key).or_default();
            entry.extend(positions);
        }
        line.clear();
    }
    for positions in map.values_mut() {
        positions.sort_unstable();
        positions.dedup();
    }
    Ok(map)
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    index: &Path,
    reads: &Path,
    read_list: &Path,
    out: Option<&Path>,
    hits_out: Option<&Path>,
    hits_intersection_out: Option<&Path>,
    gene_map: Option<&Path>,
    minimizer_positions: Option<&Path>,
    minimizer_out: Option<&Path>,
    positions_visited_out: Option<&Path>,
    local_kmer_out: Option<&Path>,
    kallisto_enum: bool,
    kallisto_strict: bool,
    kallisto_local_fallback: bool,
    kallisto_fallback: bool,
    discard_special_only: bool,
    skip_overcrowded_minimizer: bool,
    kallisto_direct_kmer: bool,
    kallisto_bifrost_find: bool,
    strand: kallistors::pseudoalign::Strand,
    fragment_length: Option<f64>,
    single_overhang: bool,
    min_range: usize,
    do_union: bool,
    no_jump: bool,
    dfk_onlist: bool,
    fr_stranded: bool,
    rf_stranded: bool,
) -> Result<()> {
    if fr_stranded && rf_stranded {
        return Err(anyhow!(
            "--fr-stranded and --rf-stranded are mutually exclusive"
        ));
    }
    let strand_specific = if fr_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::FR)
    } else if rf_stranded {
        Some(kallistors::pseudoalign::StrandSpecific::RF)
    } else {
        None
    };
    let options = kallistors::pseudoalign::PseudoalignOptions {
        min_range,
        do_union,
        dfk_onlist,
        strand_specific,
        no_jump,
        kallisto_enum,
        kallisto_strict,
        kallisto_local_fallback,
        kallisto_fallback,
        discard_special_only,
        skip_overcrowded_minimizer,
        kallisto_direct_kmer,
        kallisto_bifrost_find,
        bias: false,
        max_bias: 0,
    };

    let filter =
        fragment_length
            .filter(|v| *v > 0.0)
            .map(|v| kallistors::pseudoalign::FragmentFilter {
                fragment_length: v as u32,
                single_overhang,
            });

    let mut reader = kallistors::io::open_fastq_reader(reads)?;
    let index = if kallisto_direct_kmer {
        kallistors::pseudoalign::build_bifrost_index_with_kmer_pos(index, false)?
    } else if kallisto_fallback {
        kallistors::pseudoalign::build_bifrost_index_with_kmer(index)?
    } else {
        kallistors::pseudoalign::build_bifrost_index(index)?
    };
    let gene_map = match gene_map {
        Some(path) => Some(load_gene_map(path)?),
        None => None,
    };
    let minimizer_positions = match minimizer_positions {
        Some(path) => Some(load_minimizer_positions(path)?),
        None => None,
    };
    if minimizer_positions.is_some() && minimizer_out.is_none() {
        return Err(anyhow!("--minimizer-positions requires --minimizer-out"));
    }

    let read_list_file = File::open(read_list)?;
    let mut targets: HashSet<String> = HashSet::new();
    for line in BufReader::new(read_list_file).lines() {
        let line = line?;
        let key = normalize_header(&line);
        if !key.is_empty() {
            targets.insert(key);
        }
    }
    if targets.is_empty() {
        return Err(anyhow!("read list is empty"));
    }

    let writer: Box<dyn Write> = match out {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    let mut writer = writer;
    let mut hits_writer: Option<Box<dyn Write>> = match hits_out {
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    let mut intersection_writer: Option<Box<dyn Write>> = match hits_intersection_out {
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    let mut minimizer_writer: Option<Box<dyn Write>> = match minimizer_out {
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    let mut positions_visited_writer: Option<Box<dyn Write>> = match positions_visited_out {
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    let mut local_kmer_writer: Option<Box<dyn Write>> = match local_kmer_out {
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    if gene_map.is_some() {
        writeln!(
            writer,
            "read\taligned\treason\tec_before\tec_before_genes\tec_after_fragment\tec_after_fragment_genes\tec_after_strand\tec_after_strand_genes\tkmer_pos\tmin_pos\tkmer_seq\tmin_seq\tpositions\tsample_positions\tused_revcomp"
        )?;
    } else {
        writeln!(
            writer,
            "read\taligned\treason\tec_before\tec_after_fragment\tec_after_strand\tkmer_pos\tmin_pos\tkmer_seq\tmin_seq\tpositions\tsample_positions\tused_revcomp"
        )?;
    }
    if let Some(hits_writer) = hits_writer.as_mut() {
        if gene_map.is_some() {
            writeln!(
                hits_writer,
                "read\thit_idx\tunitig_id\tblock_idx\tread_pos\tused_revcomp\tec_size\tofflist\tec\tgene_names"
            )?;
        } else {
            writeln!(
                hits_writer,
                "read\thit_idx\tunitig_id\tblock_idx\tread_pos\tused_revcomp\tec_size\tofflist\tec"
            )?;
        }
    }
    if let Some(intersection_writer) = intersection_writer.as_mut() {
        if gene_map.is_some() {
            writeln!(
                intersection_writer,
                "read\thit_idx\tunitig_id\tblock_idx\tread_pos\tused_revcomp\tec_size\tintersection_size\tintersection_ec\tintersection_gene_names"
            )?;
        } else {
            writeln!(
                intersection_writer,
                "read\thit_idx\tunitig_id\tblock_idx\tread_pos\tused_revcomp\tec_size\tintersection_size\tintersection_ec"
            )?;
        }
    }
    if let Some(minimizer_writer) = minimizer_writer.as_mut() {
        writeln!(
            minimizer_writer,
            "read\tread_pos\tkmer\tkmer_canon\tused_revcomp\tmin_pos\tmin_seq\tmphf_key\tmphf_level\tmphf_hash_raw\tmphf_hash_domain\tmphf_bucket\tmphf_rank\tmphf_final_hash_hit\tmphf_hit\tpositions\tsample_positions\tmatched\tmatch_unitig_id\tmatch_unitig_pos\tmatch_block_idx\tis_km\tspecial_code\tspecial_in_dlist\tspecial_in_hmap\tspecial_uid\tnote"
        )?;
    }
    if let Some(writer) = positions_visited_writer.as_mut() {
        writeln!(writer, "read\tpositions_visited")?;
    }
    if let Some(local_kmer_writer) = local_kmer_writer.as_mut() {
        writeln!(
            local_kmer_writer,
            "read\tread_pos\tkmer\tlocal_hit\tlocal_ec_size"
        )?;
    }

    let mut found = 0usize;
    while let Some(record) = reader.next_record() {
        let record = record?;
        let header = String::from_utf8_lossy(&record.header).into_owned();
        let key = normalize_header(&header);
        if !targets.contains(&key) {
            continue;
        }
        if let Some(local_kmer_writer) = local_kmer_writer.as_mut() {
            let hits = kallistors::pseudoalign::local_kmer_hits(
                &index,
                &record.seq,
                kallisto_enum,
                kallisto_strict,
                skip_overcrowded_minimizer,
                kallisto_bifrost_find,
                256,
            );
            for (pos, kmer, hit, ec_size) in hits {
                writeln!(
                    local_kmer_writer,
                    "{}\t{}\t{}\t{}\t{}",
                    key,
                    pos,
                    kmer,
                    if hit { 1 } else { 0 },
                    ec_size
                )?;
            }
        }
        let (trace, hits) = kallistors::pseudoalign::trace_read_bifrost_with_hits(
            &index,
            &record.seq,
            strand,
            filter,
            options,
        );
        if let Some(writer) = positions_visited_writer.as_mut() {
            if let Some(positions) = trace.positions_visited.as_ref() {
                let joined = positions
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(writer, "{}\t{}", header, joined)?;
            } else {
                writeln!(writer, "{}\t-", header)?;
            }
        }
        let aligned = trace.reason.is_none();
        let reason = trace
            .reason
            .map(|r| r.to_string())
            .unwrap_or_else(|| "ok".to_string());
        if let Some(map) = gene_map.as_ref() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                header,
                if aligned { 1 } else { 0 },
                reason,
                format_ec(&trace.ec_before_filters),
                format_gene_names(&trace.ec_before_filters, map),
                format_ec(&trace.ec_after_fragment_filter),
                format_gene_names(&trace.ec_after_fragment_filter, map),
                format_ec(&trace.ec_after_strand_filter),
                format_gene_names(&trace.ec_after_strand_filter, map),
                trace
                    .kmer_pos
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace
                    .min_pos
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace.kmer_seq.as_deref().unwrap_or("-"),
                trace.min_seq.as_deref().unwrap_or("-"),
                trace
                    .positions
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace.sample_positions.as_deref().unwrap_or("-"),
                trace.used_revcomp
            )?;
        } else {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                header,
                if aligned { 1 } else { 0 },
                reason,
                format_ec(&trace.ec_before_filters),
                format_ec(&trace.ec_after_fragment_filter),
                format_ec(&trace.ec_after_strand_filter),
                trace
                    .kmer_pos
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace
                    .min_pos
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace.kmer_seq.as_deref().unwrap_or("-"),
                trace.min_seq.as_deref().unwrap_or("-"),
                trace
                    .positions
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                trace.sample_positions.as_deref().unwrap_or("-"),
                trace.used_revcomp
            )?;
        }
        if let Some(hits_writer) = hits_writer.as_mut() {
            for (idx, hit) in hits.iter().enumerate() {
                let ec = if hit.ec.is_empty() {
                    "-".to_string()
                } else {
                    hit.ec
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                };
                if let Some(map) = gene_map.as_ref() {
                    let genes = format_gene_names(&Some(hit.ec.clone()), map);
                    writeln!(
                        hits_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        header,
                        idx,
                        hit.unitig_id,
                        hit.block_idx,
                        hit.read_pos,
                        hit.used_revcomp,
                        hit.ec.len(),
                        hit.offlist,
                        ec,
                        genes
                    )?;
                } else {
                    writeln!(
                        hits_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        header,
                        idx,
                        hit.unitig_id,
                        hit.block_idx,
                        hit.read_pos,
                        hit.used_revcomp,
                        hit.ec.len(),
                        hit.offlist,
                        ec
                    )?;
                }
            }
        }
        if let Some(intersection_writer) = intersection_writer.as_mut() {
            let mut intersection: Option<std::collections::HashSet<u32>> = None;
            for (idx, hit) in hits.iter().enumerate() {
                let ec_set: std::collections::HashSet<u32> = hit.ec.iter().copied().collect();
                if let Some(current) = intersection.as_mut() {
                    current.retain(|v| ec_set.contains(v));
                } else {
                    intersection = Some(ec_set);
                }
                let intersection_ec = intersection
                    .as_ref()
                    .map(|set| {
                        let mut values: Vec<u32> = set.iter().copied().collect();
                        values.sort_unstable();
                        values
                    })
                    .unwrap_or_default();
                let intersection_text = if intersection_ec.is_empty() {
                    "-".to_string()
                } else {
                    intersection_ec
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                };
                if let Some(map) = gene_map.as_ref() {
                    let intersection_genes = format_gene_names(&Some(intersection_ec.clone()), map);
                    writeln!(
                        intersection_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        header,
                        idx,
                        hit.unitig_id,
                        hit.block_idx,
                        hit.read_pos,
                        hit.used_revcomp,
                        hit.ec.len(),
                        intersection_ec.len(),
                        intersection_text,
                        intersection_genes
                    )?;
                } else {
                    writeln!(
                        intersection_writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        header,
                        idx,
                        hit.unitig_id,
                        hit.block_idx,
                        hit.read_pos,
                        hit.used_revcomp,
                        hit.ec.len(),
                        intersection_ec.len(),
                        intersection_text
                    )?;
                }
            }
        }
        if let (Some(minimizer_writer), Some(minimizer_positions)) =
            (minimizer_writer.as_mut(), minimizer_positions.as_ref())
            && let Some(pos_list) = minimizer_positions.get(&key)
        {
            let rows = kallistors::pseudoalign::debug_minimizer_lookup(
                &index,
                &record.seq,
                strand,
                pos_list,
                kallisto_bifrost_find,
                kallisto_strict,
                skip_overcrowded_minimizer,
            );
            for row in rows {
                writeln!(
                    minimizer_writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    header,
                    row.read_pos,
                    row.kmer,
                    row.kmer_canon,
                    row.used_revcomp,
                    row.min_pos
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.min_seq.as_deref().unwrap_or("-"),
                    row.mphf_key.as_deref().unwrap_or("-"),
                    row.mphf_level
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.mphf_hash_raw
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.mphf_hash_domain
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.mphf_bucket
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.mphf_rank
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.mphf_final_hash_hit
                        .map(|v| if v { "1" } else { "0" })
                        .unwrap_or("-"),
                    if row.mphf_hit { 1 } else { 0 },
                    row.positions
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.sample_positions.as_deref().unwrap_or("-"),
                    if row.matched { 1 } else { 0 },
                    row.match_unitig_id
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.match_unitig_pos
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.match_block_idx
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.is_km.map(|v| if v { "1" } else { "0" }).unwrap_or("-"),
                    row.special_code
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.special_in_dlist
                        .map(|v| if v { "1" } else { "0" })
                        .unwrap_or("-"),
                    row.special_in_hmap
                        .map(|v| if v { "1" } else { "0" })
                        .unwrap_or("-"),
                    row.special_uid
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string()),
                    row.note
                )?;
            }
        }

        targets.remove(&key);
        found += 1;
        if targets.is_empty() {
            break;
        }
    }

    if !targets.is_empty() {
        eprintln!("warning: {} read(s) not found in FASTQ", targets.len());
    }
    eprintln!("traced {} read(s)", found);
    Ok(())
}
