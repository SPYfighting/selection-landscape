"""Full-run reconstruction preserves standard paired FASTQ read identities."""

import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

from selection_landscape.ngs.rebuild import count_fastq, paired_fastq_records


_TOOLS_AVAILABLE = all(shutil.which(name) for name in ("bowtie2", "bowtie2-build", "samtools"))
_LIBRARIES_AVAILABLE = all(importlib.util.find_spec(name) for name in ("pysam", "Bio"))


class PairedReadIdentityTests(unittest.TestCase):
    def test_canonical_read_name_collision_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="read-identity-") as temporary:
            first, second = Path(temporary) / "R1.fastq", Path(temporary) / "R2.fastq"
            first.write_text("@read/1\nACGT\n+\nIIII\n@read\nACGT\n+\nIIII\n")
            second.write_text("@read/2\nACGT\n+\nIIII\n@read\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "duplicate.*identifier"):
                list(paired_fastq_records(first, second))


@unittest.skipUnless(_TOOLS_AVAILABLE and _LIBRARIES_AVAILABLE,
                     "requires Bowtie2, samtools, pysam and Biopython")
class FullRunReadIdentityTests(unittest.TestCase):
    def test_mate_suffixes_preserve_full_run_molecule_counts(self):
        with tempfile.TemporaryDirectory(prefix="read-identity-") as temporary:
            root = Path(temporary)
            # Two synthetic UID families, four paired reads per family.
            codons = ["GAA", "GTT", "CAA", "GCT", "TGC", "GAC", "TTC", "GGT",
                      "CAT", "ATT", "AAA", "CTG", "ATG", "AAC", "CCT", "CGT",
                      "TCT", "ACT", "TGG", "TAT", "GAC", "GGT", "TTC", "AAC",
                      "ATG", "CAT", "CCT", "AAA", "CGT", "GCT", "TCT", "CTG",
                      "TGG", "CAA", "TGC", "ATT", "GTT", "ACT", "GTT", "TCT", "TCT"]
            coding = "".join(codons)
            anchor = "ATAAAACGACGGC"
            reference = root / "reference.fa"
            reference.write_text(">synthetic\n" + coding + "\n", encoding="ascii")
            configuration = {
                "uid": {"anchor": anchor, "uid_length": 14},
                "orf": {"nt_start_zero_based": 0, "nt_end_exclusive": len(coding),
                        "frame_zero_based": 0, "aa_length": len(codons),
                        "start_motif": "EVQ", "end_motif": "VTVSS", "lead_aa": "unused"},
                "parameter_sets": [{"name": "strict", "anchor_mismatch_max": 0,
                                    "bowtie2_mode": "--local --very-sensitive-local",
                                    "min_mapq": 10, "min_orf_coverage": 1.0}],
            }
            results = []
            for suffix in (False, True):
                style = root / ("suffixed" if suffix else "plain")
                style.mkdir()
                first, second = style / "R1.fastq", style / "R2.fastq"
                with first.open("w") as r1, second.open("w") as r2:
                    for index in range(8):
                        uid = ("ACGATCGTACGTAC", "TGCTAGCATGCTAG")[index // 4]
                        sequence = coding + anchor + uid
                        reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                        for handle, mate, bases in ((r1, 1, sequence), (r2, 2, reverse)):
                            identifier = f"read{index}" + (f"/{mate}" if suffix else "")
                            handle.write(f"@{identifier}\n{bases}\n+\n{'I' * len(bases)}\n")
                manifest = style / "samples.csv"
                with manifest.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["sample_id", "round", "chain", "r1_path", "r2_path"])
                    writer.writeheader()
                    writer.writerow(dict(sample_id="synthetic", round=1, chain="Crick",
                                         r1_path=first.name, r2_path=second.name))
                summary = count_fastq(manifest, reference, configuration, style / "output",
                                      parameter_name="strict", max_reads=None, threads=1)[0]
                results.append({key: summary[key] for key in
                                ("reads_in", "uid_found", "aligned", "molecule_count")})
            expected = {"reads_in": 8, "uid_found": 8, "aligned": 16, "molecule_count": 2}
            self.assertEqual(results, [expected, expected])


if __name__ == "__main__":
    unittest.main()
