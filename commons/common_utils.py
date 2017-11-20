'''
borrow from BioTools
clean variant id
'''

'''
request ensembl for bases based on location
'''
def find_bases(chrom,start,end=None,build='hg19',strand=1):
    import requests
    import time
    # translate build
    if build=='hg19':
        server = "http://grch37.rest.ensembl.org"
    elif build=='hg38':
        server = "http://rest.ensembl.org"
    end = end or start
    ext = '''/sequence/region/human/%(chrom)s:%(start)s..%(end)s:%(strand)s''' % locals()
    attempt = 5
    while attempt:
        try:
            r = requests.get(server+ext, headers={'Content-Type':'application/json' })
            time.sleep(0.05)
            if r.ok:
                break
        except requests.HTTPError:
            print('query ensembl HTTPError, retry')
            attempt -= 1
            time.sleep(2)
        except requests.ConnectionError:
            print('query ensembl ConnectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    '''
    if not r.ok:
        return r.raise_for_status()
    '''
    decoded = r.json()
    return str(decoded['seq'])

def clean_variant(v,fasta_ref=None,build='hg19'):
    # sometimes variant has funny format, which has more - than expected, such as 1-117122294---TCT.
    # since currently I can't find where fasta_ref is,
    #  use find_bases to fill in the gap
    #  until I can use samtools to query fasta_ref
    if v.count('-') == 4:
        import subprocess
        if v[-1] == '-':
            # deletion
            chrom,pos,ref,rubbish,rubbish = v.split('-')
            pos = int(pos)-1
            common_base = find_bases(chrom,pos,build=build)
            ref = common_base + ref
            alt = common_base
        else:
            # insertion
            chrom,pos,ref,rubbish,alt = v.split('-')
            pos = int(pos)-1
            common_base = find_bases(chrom,pos,build=build)
            ref = common_base
            alt = common_base + alt
    else:
        chrom,pos,ref,alt = v.split('-')
        pos = int(pos)
    if len(ref) < len(alt):
        ran = range(len(ref))
    else:
        ran = range(len(alt))
    # insert
    for e in ran:
        ref_e = len(ref) - e - 1
        alt_e = len(alt) - e - 1
        if ref[ref_e] != alt[alt_e]: break
    for b in ran:
        if ref[b] != alt[b] or len(ref[b:ref_e+1]) == 1 or len(alt[b:alt_e+1]) == 1:
            break
    return '-'.join([chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]])


