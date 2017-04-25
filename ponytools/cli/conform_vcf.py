import argparse

def conform(args):
    from ponytools.MNEcAnnot import MNEc2MAnnot
    from ponytools.Exceptions import TriAllelicError
    import ponytools as pc
    # get a mapping of prones to chrom positionss
    annot = MNEc2MAnnot()
    # build the BIEC index
    id_map = {
        **{x:y for x,y in zip(annot._annot.ProbeSetID,zip(annot._annot.chrom,annot._annot.pos,annot._annot.MNEcID))},
        **{x:y for x,y in zip(annot._annot.AffySNPID,zip(annot._annot.chrom,annot._annot.pos,annot._annot.MNEcID))},
        **{x:y for x,y in zip([x.split('.')[4] for x in annot._annot.MNEcID],zip(annot._annot.chrom,annot._annot.pos,annot._annot.MNEcID))},
        **{x:y for x,y in zip(zip(annot._annot.chrom,annot._annot.pos),zip(annot._annot.chrom,annot._annot.pos,annot._annot.MNEcID))}
    }
    REF_map = {
            x:(y,z) for x,(y,z) in zip(annot._annot.MNEcID,zip(annot._annot.REF,annot._annot.ALT))
    }
     
    # Iterate over the VCF
    vcf = pc.VCF(args.vcf)
    with open(args.out,'w') as OUT:
        #Clear the header, add the actual chromosomes
        vcf.header['contig'] = []
        for chrom in annot._annot.chrom.unique():
            vcf.header['contig'].append('<ID={}>'.format(chrom))
        print(vcf.header,file=OUT)
        for i,variant in enumerate(vcf):
            # Update the chrom and position to be reflect the MNEc naming format 
            if variant.id in id_map:
                chrom,pos,mnecid = id_map[variant.id]
                variant.chrom = chrom
                variant.pos = pos
                variant.id = mnecid 
            elif ('{}'.format(variant.chrom),variant.pos) in id_map:
                chrom = variant.chrom
                pos = variant.pos
                if not chrom.startswith('chr'):
                    chrom = 'chr{}'.format(chrom)
                chrom,pos,mnecid = id_map[('{}'.format(chrom),variant.pos)]
                variant.chrom = chrom
                variant.pos = pos
                variant.id = mnecid 
            else:
                continue
            # Make sure the genotypes are conformed
            if variant.alt == '.':
                variant.alt = REF_map[variant.id][1]
            try:
                variant.conform(REF_map[variant.id][0])
            except TriAllelicError as e:
                continue
            try:
                assert variant.alt == REF_map[variant.id][1]
            except AssertionError as e:
                continue
            # Fix the borked dot genotpyes
            for i,x in enumerate(variant.genos):
                if x == '.':
                    variant.genos[i] = './.'
            print(variant,file=OUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Conform input VCFs to match MNEc Reference Population')
    parser.add_argument('--vcf')
    parser.add_argument('--fasta')
    parser.add_argument('--out')

    args = parser.parse_args()
    conform(args)

