"""Write a synthetic GERP bigWig (random piecewise-constant scores) for one contig, for the local test.
Both pipelines read the same file, so only its format matters; the values are spread so that LOFTEE's
END_TRUNC (GERP_DIST <= -58) fires for some variants and not others.

  python make_fake_gerp.py CONTIG LENGTH OUT.bw      (needs pyBigWig)
"""
import random
import sys

import pyBigWig

contig, length, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rng = random.Random(11)
starts, ends, values = [], [], []
pos = 0
while pos < length:
    end = min(length, pos + rng.randint(5, 200))
    starts.append(pos); ends.append(end)
    values.append(round(max(-12.3, min(6.2, rng.gauss(-0.8, 3.0))), 3))
    pos = end
bw = pyBigWig.open(out, "w")
bw.addHeader([(contig, length)], maxZooms=10)
bw.addEntries([contig] * len(starts), starts, ends=ends, values=values)
bw.close()
print(f"{len(starts)} intervals on {contig}", file=sys.stderr)
