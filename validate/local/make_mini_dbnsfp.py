"""Build a dbNSFP-format file for chr21 test variants, with planted edge cases for the plugin emulation.
stdin: split-vep lines 'CHROM POS REF ALT Consequence Amino_acids'. stdout: dbNSFP-like TSV (unsorted body)."""
import random
import sys

rng = random.Random(7)
INCL = {"missense_variant", "stop_lost", "stop_gained", "start_lost"}
aa_by_var = {}
for line in sys.stdin:
    chrom, pos, ref, alt, csq, aa = line.split()
    if len(ref) != 1 or len(alt) != 1 or aa == "." or "/" not in aa:
        continue
    if not INCL & set(csq.split("&")):
        continue
    aa_by_var.setdefault((chrom.replace("chr", ""), int(pos), ref, alt), [])
    if aa not in aa_by_var[(chrom.replace("chr", ""), int(pos), ref, alt)]:
        aa_by_var[(chrom.replace("chr", ""), int(pos), ref, alt)].append(aa)


def revel():
    r = rng.random()
    if r < 0.15:
        return "."
    n = rng.randint(1, 4)
    return ";".join(rng.choice([".", f"{rng.random():.3f}"]) if rng.random() < 0.3 else f"{rng.random():.3f}"
                    for _ in range(n))


def cadd():
    return "." if rng.random() < 0.1 else f"{rng.uniform(0, 45):.3f}"


rows = []
counts = dict(dup=0, dotfirst=0, decoy=0, badref=0)
for (c, p, ref, alt), aas in aa_by_var.items():
    for aa in aas:
        aref, aalt = (x.replace("*", "X") for x in aa.split("/"))
        r0 = ref
        if rng.random() < 0.05:
            r0 = rng.choice([b for b in "ACGT" if b != ref]); counts["badref"] += 1
        u = rng.random()
        if u < 0.08:  # first row both '.', a later row with values: plugin must emit nothing
            rows.append((c, p, r0, alt, aref, aalt, ".", ".")); counts["dotfirst"] += 1
            rows.append((c, p, r0, alt, aref, aalt, revel(), cadd()))
        else:
            rows.append((c, p, r0, alt, aref, aalt, revel(), cadd()))
            if u < 0.18:  # later duplicate with different values: plugin must take the first
                rows.append((c, p, r0, alt, aref, aalt, revel(), cadd())); counts["dup"] += 1
        if rng.random() < 0.1:  # decoy amino-acid change: must not match
            rows.append((c, p, ref, alt, aref, "C" if aalt != "C" else "G", revel(), cadd())); counts["decoy"] += 1

rows.sort(key=lambda r: r[1])  # stable: keeps planted within-position order
print("#chr\tpos(1-based)\tref\talt\taaref\taaalt\tREVEL_score\tCADD_phred")
for r in rows:
    print("\t".join(map(str, r)))
print(f"{len(aa_by_var)} variants, {len(rows)} rows, planted: {counts}", file=sys.stderr)
