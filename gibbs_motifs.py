from setup import *
from gibbs_sampling import run_GibbsSampler

!wget -c https://www.dropbox.com/sh/u3hp274o6cg98qn/AAA31BUmilIxZU6-bwUtQnRTa?dl=1 -O gibbs_motifs.zip
!unzip -o gibbs_motifs.zip -d gibbs_motifs
os.unlink("gibbs_motifs.zip")
gibbs_motifs = {}
for fname in os.listdir("gibbs_motifs"):
    with open(os.path.join("gibbs_motifs", fname), "r") as R:
        gibbs_motifs[fname] = R.read().strip().split()

# Test Cases, can change data
run_GibbsSampler(gibbs_motifs["data4.txt"], 10, 1)