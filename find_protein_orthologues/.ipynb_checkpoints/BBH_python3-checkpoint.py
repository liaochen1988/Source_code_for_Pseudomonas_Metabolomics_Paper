#from IPython import embed
import pandas as pd


# will DL scripts, run BBH, store etc


def dl_genome(id, folder='genomes'): # be sure get CORRECT id... 
    from glob import glob
    files =glob('%s/*.gb'%folder)
    out_file = '%s/%s.gb'%(folder, id)
    # embed()
    if out_file in files:
        print(out_file, 'already downloaded')
        #print out_file, 'already downloaded'
        return
    else:
        print('downloading %s from NCBI'%id)
        #print 'downloading %s from NCBI'%id
    from Bio import Entrez
    Entrez.email = "jmonk@ucsd.edu"     # Always tell NCBI who you are
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    fout = open(out_file,'w')
    fout.write(handle.read())
    fout.close()
    
def get_strain_info(folder='genomes'):
    from glob import glob
    from Bio import Entrez, SeqIO
    files = glob('%s/*.gb'%folder)
    strain_info = []
    for file in files:
        handle = open(file)
        print('parsing', file)
        #print 'parsing', file
        record = SeqIO.read(handle, "genbank")
        
        for f in record.features:
            if f.type=='source':
                info = {}
                info['file'] = file
                info['id'] = file.split('\\')[-1].split('.')[0]
                # embed()
                for q in f.qualifiers.keys():
                    info[q] = '|'.join(f.qualifiers[q])
                strain_info.append(info)
    return pd.DataFrame(strain_info)
            
    
def gen_gene_table(id, type='prot', in_folder='genomes', out_folder='prots'):
    from Bio import Entrez, SeqIO
    if type=='nucl':
        out_folder='nucls'
    
    in_file = '%s/%s.gbk'%(in_folder, id)
    out_file='%s/%s.csv'%(out_folder, id)
    
    from glob import glob
    files =glob('%s/*.csv'%out_folder)
    # embed()
    if out_file in files:
        print(out_file, 'table already created')
        #print out_file, 'table already created'
        return
    else:
        print('building %s table for %s'%(type, id))
        #print 'building %s table for %s'%(type, id)
        
    handle = open(in_file)
    print('parsing', in_file)
    #print 'parsing', in_file
    # record = SeqIO.read(handle, "genbank")
    records = SeqIO.parse(handle, "genbank")
    # embed()
    x = 0
    out = []
    for record in records:
        for f in record.features:
            # embed()
            if type == 'prot':
                if f.type=='CDS':
                    info = {}
                    info['type'] = f.type
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                    out.append(info)  
            elif type=='nucl':
                if f.type=='gene':
                    info = {}
                    info['type'] = f.type
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                    try:
                        info['seq'] = f.extract(record).seq.tostring()
                        info['strand'] = f.strand
                    except:
                        pass
                    out.append(info) 
    out=pd.DataFrame(out)
    out.to_csv(out_file)

def parse_genome(id, type='prot', in_folder='genomes', out_folder='prots'):
    from Bio import Entrez, SeqIO
    if type=='nucl':
        out_folder='nucls'
    
    in_file = '%s/%s.gbk'%(in_folder, id)
    out_file='%s/%s.faa'%(out_folder, id)
    
    from glob import glob
    files =glob('%s/*.faa'%out_folder)
    if out_file in files:
        print(out_file, 'already parsed')
        #print out_file, 'already parsed'
        return
    else:
        print('parsing %s'%id)
        #print 'parsing %s'%id
    
    handle = open(in_file)
    print('parsing', in_file)
    #print 'parsing', in_file
    
    fout = open(out_file,'w')
    x = 0
    
    records = SeqIO.parse(handle, "genbank")
    # record = SeqIO.read(handle, "genbank")
    for record in records:
        # print record
        for f in record.features:
            if f.type=='CDS':
                seq=f.extract(record.seq)
                # embed()
                if type=='nucl':
                    seq=seq.tostring()
                else:
                    seq=seq.translate().tostring()
                if 'locus_tag' in f.qualifiers.keys():
                    locus = f.qualifiers['locus_tag'][0]
                elif 'gene' in f.qualifiers.keys():
                    locus = f.qualifiers['gene'][0]
                else:
                    locus = x
                    x+=1
                # print locus, seq
                fout.write('>%s\n%s\n'%(locus, seq))
    fout.close()
    
    


def make_blast_db(id,db_type='prot', overwrite=False):
    import os
    if os.path.isfile(id+'.phr') and overwrite==False:
        print(id, db_type,' blast db exists, skipping')
        #print id, db_type,' blast db exists, skipping'
        return
    # out_file ='%s\\%s.fa.pin'%(folder, id)
    # from glob import glob
    # files =glob('%s/*.fa.pin'%folder)
    # if out_file in files:
        # print id, 'already has a blast db'
        # return
    if db_type=='prot':
        infile = 'prots/%s.faa'%(id)
    elif db_type=='nucl':
        infile = 'nucls/%s.faa'%(id)
    
    """ make a blast db (BLAST 2.2.26+) """
    cmd_line='makeblastdb -in %s -dbtype %s' %(infile,db_type)
    print('making blast db with following command line...')
    print(cmd_line)
    os.system(cmd_line)

def run_blastp(seq,db,in_folder='prots', out_folder='bbh', \
    out=None,outfmt=6,evalue=0.001,threads=12):
    import os
    if out==None:
        out='%s/%s_vs_%s.txt'%(out_folder, seq, db)
        
    
    from glob import glob
    files =glob('%s/*.txt'%out_folder)
    if out in files:
        print(seq, 'already blasted')
        return
    
    print('blasting %s vs %s'%(seq, db))
    
    db = '%s/%s.faa'%(in_folder, db)
    seq = '%s/%s.faa'%(in_folder, seq)
    cmd_line='blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print('running blastp with following command line...')
    print(cmd_line)
    os.system(cmd_line)
    return out
    
def run_blastn(seq,db,out=None,outfmt=6,evalue=0.001,threads=1):
    import os
    # if out==None:
        # out='%s\\%s_vs_%s_nucl.txt'%(out_folder, seq, db)
        
    # from glob import glob
    # files =glob('%s/*.txt'%out_folder)
    # if out in files:
        # print seq, 'already blasted'
        # return
    
    print('blasting %s vs %s'%(seq, db))
    
    db = '%s'%(db)
    seq = '%s'%(seq)
    cmd_line='blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print('running blastn with following command line...')
    print(cmd_line)
    os.system(cmd_line)
    return out
    
def get_gene_lens(query, in_folder='prots'):
    from Bio import SeqIO
    file = '%s/%s.faa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out
    

def get_bbh(query, subject, in_folder='bbh', outfile=''):
    # dl_genome(query)
    # parse_genome(query)
    # gen_gene_table(id, type='prot')
    # gen_gene_table(id, type='nucl')
    # make_blast_db(query)
    
    # dl_genome(subject)
    # parse_genome(subject)
    # gen_gene_table(id, type='prot')
    # gen_gene_table(id, type='nucl')
    # make_blast_db(subject)
    
    run_blastp(query, subject)
    run_blastp(subject, query)
    
    query_lengths = get_gene_lens(query, in_folder='prots')
    subject_lengths = get_gene_lens(subject, in_folder='prots')
    
    # run_blastn(query, subject)
    # run_blastn(subject, query)
    if outfile=='':
        outfile = '%s/%s_vs_%s_parsed.csv'%(in_folder,query, subject)
    
    from glob import glob
    files=glob('%s/*_parsed.csv'%in_folder)
    
    # if outfile in files:
        # print 'bbh already parsed for', query, subject
        # out = pd.read_csv(outfile)
        # return out
    print('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,query, subject), sep='\t', names=cols)
    bbh = pd.merge(bbh, query_lengths) 
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']
    bbh2=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,subject, query), sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths) 
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()
    # embed()
    
    # FILTER GENES THAT HAVE < 80% coverage
    bbh = bbh[bbh.COV>=0.8]
    bbh2 = bbh2[bbh2.COV>=0.8]
    
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.loc[res.PID.idxmax()]
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        # embed()
        if g==best_gene2:
            # print g, '<=>', best_gene
            best_hit['BBH'] = '<=>'
        else:
            # print g, '->', best_gene
            best_hit['BBH'] = '->'
            # embed()
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])
    # embed()
    out.to_csv(outfile)
    return out
    
def convert_to_fasta(file):
    outfile = file[:-3]+'fasta'
    fout = open(outfile,'w')
    from Bio import SeqIO
    handle=open(file)
    out = []
    Ncount = 0
    for record in SeqIO.parse(handle,"genbank"):
        # print record.seq
        fout.write('>%s\n%s\n'%(record.id, record.seq))
        for s in record.seq:
            out.append(s)
            if s=='N':
                Ncount+=1
        out = list(set(out))
    
    fout.close()
    return outfile
    
def get_MLST(strain, mlst):
    import os
    ''' takes a gbk file, converts to fasta, creates a blast dB and blasts the fasta-formatted mlst file against it '''
    
    fasta = strain[:-3]+'fasta'
    if not os.path.isfile(fasta):
        print('converting %s to fasta'%strain)
        fasta = convert_to_fasta(strain)

    if not os.path.isfile(fasta+'.nhr'):
        print('making blast db')
        make_blast_db(fasta,db_type='nucl')
    
    run_blastn(mlst,fasta, out='mlst_tmp.csv')
    
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('mlst_tmp.csv', sep='\t', names=cols)
    if len(bbh)>0:
        max_bit = bbh.bitScore.max()
        max_bit_idx = bbh.bitScore.idxmax()
        # print bbh.loc[max_bit_idx]
        
        return bbh.loc[max_bit_idx]['gene']
    else:
        return 'error'
        
        
def get_MLST(fasta, mlst):
    import os
    ''' takes a gbk file, converts to fasta, creates a blast dB and blasts the fasta-formatted mlst file against it '''
    
    # fasta = strain[:-3]+'fasta'
    if not os.path.isfile(fasta):
        print('converting %s to fasta'%fasta)
        fasta = convert_to_fasta(fasta)

    if not os.path.isfile(fasta+'.nhr'):
        print('making blast db')
        make_blast_db(fasta,db_type='nucl')
    
    run_blastn(mlst,fasta, out='mlst_tmp.csv')
    
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('mlst_tmp.csv', sep='\t', names=cols)
    if len(bbh)>0:
        max_bit = bbh.bitScore.max()
        max_bit_idx = bbh.bitScore.idxmax()
        print(bbh.loc[max_bit_idx])
        # embed()
        
        return bbh.loc[max_bit_idx]['gene']
    else:
        return 'error'
    
    
if __name__=='__main__':

    # helping pekar:    
    # file = 'C_ill_clc_rast_cfb.gbk'
    # 
    
    # getting best blast hit from MLST seqs to NA of each genome
    errors = []
    # 1 build dBs
    from glob import glob
    # strain = '../genomes/CD_BTH-0054.gbf'
    strains = glob('../genomes/*.gbf')
    mlst_files = glob('MLST_profiles/*.fas')
    all_out = {}
    for strain in strains:
        out = {}
        strain_name = strain.split('\\')[-1].split('.')[0]
        print(strain_name)
        for mlst_file in mlst_files:
            res = get_MLST(strain, mlst_file)
            if len(res.split('_'))>1:
                gene = res.split('_')[0]
                out[gene]=res.split('_')[1]
            else:
                out['error']='error'
        all_out[strain_name]=out
        print(strain_name, out)
        # embed()
    data = pd.DataFrame(all_out).transpose()
    del data['error']
    data = data.dropna().reset_index().rename(columns={'index':'strain'})
    data.to_csv('mlst_results.csv')
    
    to_compare = 'adk  atpA  dxr  glyA  recA  sodA  tpi'
    to_compare = to_compare.split('  ')
    for c in to_compare:
        data[c] = pd.to_numeric(data[c])
        # to_compare.append(c)
    
    # parse mlst results
    
    mlsts = pd.read_csv('bigsdb.txt',sep='\t')
    
    res = pd.merge(data, mlsts,left_on=to_compare, right_on=to_compare,how='left')
        
        
    #embed()
    
    

