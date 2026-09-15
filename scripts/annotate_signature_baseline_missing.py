"""Complete baseline-only annotation without requiring a post-treatment sample."""
from pathlib import Path
import annotate_cd8_reference as base

def main():
    base.OUT = base.ROOT / 'results/signature_revision_20260915/additional_annotations'
    base.OBJECTS = base.ROOT / 'data/processed/signature_revision_20260915'
    base.OUT.mkdir(parents=True, exist_ok=True)
    base.OBJECTS.mkdir(parents=True, exist_ok=True)
    low=base.Model.load(str(base.ROOT/'models/celltypist/Immune_All_Low.pkl'))
    high=base.Model.load(str(base.ROOT/'models/celltypist/Immune_All_High.pkl'))
    with base.threadpool_limits(limits=4):
        for patient in ['HN17','HN31']:
            base.annotate_sample(patient,'pre',low,high)

if __name__=='__main__': main()
