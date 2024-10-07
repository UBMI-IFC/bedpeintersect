#!/usr/bin/env python3

import argparse
import re
from sys import stdout
from pathlib import Path
import pandas as pd
import numpy as np
from pybedtools import BedTool

genename =  re.compile(r'gene_name "(?P<gene>.*?)"')
biotype = re.compile(r'gene_biotype "(?P<biotipye>.*?)"')
geneid = re.compile(r'gene_id "(?P<geneid>.*?)"')
exonnumber = re.compile(r'exon_number "(?P<exonnumber>.*?)"')

def argumetns():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedpefile',
                        help='HiC bedpe contact file')
    parser.add_argument('bedsfolder',
                        help='Folder with bed files to intersect')
    return parser.parse_args()

def melt_tag_beddpee(bedpefile, outfile=None):
    """Transform a bedpee file into a melted version. Where each pair of a
    contact is stored in a different line"""
    myfile = Path(bedpefile)
    # Outfile
    if outfile:
        output = open(outfile, 'w')
    else:
        output = stdout
    # count = 1
    countname = 1
    with open(myfile) as inf, output:
        for line in inf:
            if line[0] == '#' or line == '':
                continue
            contact = f'cont_{countname:0>5}'
            parts = line.split()
            a = '\t'.join(parts[:3]) + f'\t{contact}a\n'
            b = '\t'.join(parts[3:]) + f'\t{contact}b\n'
            output.write(a)      # a contact
            output.write(b)      # b contact
            # add to counter
            countname += 1



if __name__ == '__main__':
    # input files
    bedpefile = 'P046T_neoloops_095_mres.bed'
    bedpemeltfile = 'melt_' + bedpefile
    bedspath = Path('./beds')
    bedfiles = list(bedspath.rglob('*.bed'))
    beds = [BedTool(f) for f in bedfiles]

    # metl bedpe
    melt_tag_beddpee(bedpefile, bedpemeltfile)
    melted = BedTool(bedpemeltfile)

    # Intersect
    intersect = melted.intersect([bed.fn for bed in beds],
                                 filenames=True,
                                 wa=True,
                                 C=True)

    # Create df
    df_inter = pd.read_table(intersect.fn,
                             names=['chr', 'init', 'end', 'id', 'filename',
                                    'counts'],
                             index_col=3)
    df = df_inter.pivot(columns=['filename'], values='counts')
    ranges = df_inter[['chr', 'init', 'end']].iloc[::2]
    df = ranges.join(df)

    # closest gene
    melt_sorted = melted.sort()
    closegene = melt_sorted.closest('genes_chr_sorted.gtf.gz', d=True, t='first')
    closeexon = melt_sorted.closest('exon_chr_sorted.gtf.gz', d=True, t='first')
    closeexoncomplete = melt_sorted.closest('exon_chr_sorted.gtf.gz', d=True, t='all')
    assert len(closegene) == len(closeexon), 'Annotation intersection wrong sizes'


    df['gene name'] = np.nan
    df['gene id'] = np.nan
    df['gene biotype'] = np.nan
    df['dist'] = np.nan
    df['region'] = np.nan

    condits = [0,0,0]
    for i, (gene, part) in enumerate(zip(closegene, closeexon)):
        fields = gene.fields
        partfields = part.fields
        attr = fields[-2]
        attrpart = partfields[-2]
        distance = int(fields[-1])
        partdistance = int(partfields[-1])
        gn = genename.findall(attr)
        gnp = genename.findall(attrpart)
        gi = geneid.findall(attr)
        gip = geneid.findall(attrpart)
        parttype = attrpart[6]
        gb = biotype.findall(attr)
        gbp = biotype.findall(attrpart)

        df.loc[gene.name, 'gene name'] = gn[0]
        df.loc[gene.name, 'gene id'] = gi[0]
        df.loc[gene.name, 'gene biotype'] = gb[0]
        df.loc[gene.name, 'dist'] = distance

        # The following code for gene region specification
        # is likely inestable. Mean reasons
        #  1. HiC regions are large, many features can fint
        #  2. Mixed annotations for a region require many rules to
        #    accurately map it
        if distance > 0:
            # not in a gene
            df.loc[gene.name, 'region'] = 'intergenic'
            condits[0] += 1
        elif gi[0] == gip[0]:
            condits[1] += 1
            # in a gene and easy exon recognicing
            if 'utr' in partfields[6]:
                df.loc[gene.name, 'region'] = partfields[6]
            elif partfields[6] == 'exon':
                # n = exonnumber.findall(attrpart)[0]
                df.loc[gene.name, 'region'] = 'exon' # + n
            else:
                break
        else:
            condits[2] += 1
            genf = BedTool([gene])
            features = closeexoncomplete.intersect(genf)
            genes = set([geneid.findall(a.fields[-2])[0] for a in features])
            if gi[0] in genes:
                df.loc[gene.name, 'region'] = 'exon'
            else:
                df.loc[gene.name, 'region'] = 'intron'

    # Write bedpee file
    with open('pru.bedpe', 'w') as outbedpe:
        for i in range(len(df)):
            if i % 2 == 0:
                contcta = [str(j) for j in df.iloc[i]]
            else:
                contctb = [str(j) for j in df.iloc[i]]
                contactname = df.index[i][:-1]
                # Line construction
                line = '\t'.join(contcta[:3]) + '\t'
                line += '\t'.join(contctb[:3]) + '\t'
                line += contactname + '\t'
                lenbeds = len(beds)
                line += '\t'.join(contcta[3:3+lenbeds]) + '\t'
                line += '\t'.join(contctb[3:3+lenbeds]) + '\n'
                outbedpe.write(line)
