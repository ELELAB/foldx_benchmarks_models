'''

Created on Apr 17, 2018
@author: Tycho

Calculate the sequence identities from a hhpred .hhr alignment file.
'''

import os
import urllib2
from Bio import SeqIO
from Bio.PDB import *
from StringIO import StringIO
from operator import itemgetter
import re
import time

def main():
        
#     os.chdir("Z:\mutatex_project\P62942")
#     ali_file = "hhpred_P62942.hhr"
#     query_fasta = "P62942.fasta"

    os.chdir("Z:\mutatex_project\P0CG48")
    ali_file = "hhpred_Ubiquetin.hhr"
    query_fasta = "P0CG48.fasta"

#     os.chdir("Z:\mutatex_project\P61626")
#     ali_file = "hhpred_P61626.hhr"
#     query_fasta = "P61626.fasta"

#     os.chdir("Z:\mutatex_project\FKB1A")
#     ali_file = "hhpred_FKB1A.hhr"
#     query_fasta = "FKB1A.fasta"

    
    with open(query_fasta, 'r') as f:
        raw_query_seq = ('')
        for line in f:
            if not line.startswith(">") and line is not "":
                raw_query_seq += line.replace("\n", "")

    seqs = get_sequences(ali_file)

    ninety_seventy = {}
    seventy_fifty = {}
    fifty_thirty = {}
    thirty_twenty = {}
    for seq in seqs:
        alignment = ["." if x is not y else "|" for x, y in zip(seq[1][0], seq[1][1]) ]
        seq_identity = len(raw_query_seq)
        print seq[1][0]
        for x, y in zip(seq[1][0], seq[1][1]):
            if x is not y:
                seq_identity -= 1
        print ""
        seq_id_percent = float(seq_identity) / len(raw_query_seq) * 100
        if seq_identity != 0 and len(seq[1][1]) > float(len(raw_query_seq)) * 0.85:
            print seq[0] 
            print "Len query seq: " + str(len(raw_query_seq))
            print "Len template seq " + str(len(seq[1][1]))
            print "Seq identity: " + str(seq_identity)
            print "Seq identity: " + format(seq_id_percent, '.2f') + "%"
            print seq[1][0]
            print "".join(alignment)
            print seq[1][1]
            print ""
            if(seq_id_percent < 90 and seq_id_percent >= 70):
                ninety_seventy[seq[0].split("_")[0]] = {"group" : "90-70",
                                                    "pdb_id": seq[0].split("_")[0],
                                                    "chain" : seq[0].split("_")[-1],
                                                    "query_len" : len(raw_query_seq),
                                                    "template_len" : len(seq[1][1]),
                                                    "seq_identity" : seq_identity,
                                                    "seq_identity_percent" : format(seq_id_percent, '.2f'),
                                                    "template_ali_start" : seq[2][0],
                                                    "template_ali_end" : seq[2][1],
                                                    "query_sequence" : seq[1][0],
                                                    "template_sequence" : seq[1][1],
                                                    "alignment" : "".join(alignment),
                                                    }
            elif(seq_id_percent < 70 and seq_id_percent >= 50):
                seventy_fifty[seq[0].split("_")[0]] = {"group" : "70-50",
                                                    "pdb_id": seq[0].split("_")[0],
                                                    "chain" : seq[0].split("_")[-1],
                                                    "query_len" : len(raw_query_seq),
                                                    "template_len" : len(seq[1][1]),
                                                    "seq_identity" : seq_identity,
                                                    "seq_identity_percent" : format(seq_id_percent, '2f'),
                                                    "template_ali_start" : seq[2][0],
                                                    "template_ali_end" : seq[2][1],
                                                    "query_sequence" : seq[1][0],
                                                    "template_sequence" : seq[1][1],
                                                    "alignment" : "".join(alignment),
                                                    }
            elif(seq_id_percent < 50 and seq_id_percent >= 30):
                fifty_thirty[seq[0].split("_")[0]] = {"group" : "50-30",
                                                    "pdb_id": seq[0].split("_")[0],
                                                    "chain" : seq[0].split("_")[-1],
                                                    "query_len" : len(raw_query_seq),
                                                    "template_len" : len(seq[1][1]),
                                                    "seq_identity" : seq_identity,
                                                    "seq_identity_percent" : format(seq_id_percent, '2f'),
                                                    "template_ali_start" : seq[2][0],
                                                    "template_ali_end" : seq[2][1],
                                                    "query_sequence" : seq[1][0],
                                                    "template_sequence" : seq[1][1],
                                                    "alignment" : "".join(alignment),
                                                    }
            elif(seq_id_percent < 30 and seq_id_percent >= 20):
                thirty_twenty[seq[0].split("_")[0]] = {"group" : "30-20",
                                                    "pdb_id": seq[0].split("_")[0],
                                                    "chain" : seq[0].split("_")[-1],
                                                    "query_len" : len(raw_query_seq),
                                                    "template_len" : len(seq[1][1]),
                                                    "seq_identity" : seq_identity,
                                                    "seq_identity_percent" : format(seq_id_percent, '2f'),
                                                    "template_ali_start" : seq[2][0],
                                                    "template_ali_end" : seq[2][1],
                                                    "query_sequence" : seq[1][0],
                                                    "template_sequence" : seq[1][1],
                                                    "alignment" : "".join(alignment),
                                                    }
        else:
            print seq[0] + "\nSeq identity: 0"
            
    
    print len(ninety_seventy)
    print ninety_seventy.keys()
    print len(seventy_fifty)
    print seventy_fifty.keys()
    print len(fifty_thirty)
    print fifty_thirty.keys()
    print len(thirty_twenty)
    print thirty_twenty.keys()
    
    pdb_downloader(ninety_seventy, "90-70")
    pdb_downloader(seventy_fifty, "70-50")
    pdb_downloader(fifty_thirty, "50-30")
    pdb_downloader(thirty_twenty, "30-20")
       
    pdb_parser(ninety_seventy, "90-70")
    pdb_parser(seventy_fifty, "70-50")     
    pdb_parser(fifty_thirty, "50-30")
    pdb_parser(thirty_twenty, "30-20")
     
    ninety_seventy_list = []
    seventy_fifty_list = []
    fifty_thirty_list = []
    thirty_twenty_list = []
     
    for key, value in ninety_seventy.iteritems():
        if 'resolution' not in value:
            value['resolution'] = 0
        ninety_seventy_list.append(value)
    for key, value in seventy_fifty.iteritems():
        if 'resolution' not in value:
            value['resolution'] = 0
        seventy_fifty_list.append(value)
    for key, value in fifty_thirty.iteritems():
        if 'resolution' not in value:
            value['resolution'] = 0
        fifty_thirty_list.append(value)
    for key, value in thirty_twenty.iteritems():
        if 'resolution' not in value:
            value['resolution'] = 0
        thirty_twenty_list.append(value)
     
    ninety_seventy_list  = sorted(ninety_seventy_list, key=itemgetter('resolution'))
    seventy_fifty_list = sorted(seventy_fifty_list, key=itemgetter('resolution'))
    fifty_thirty_list = sorted(fifty_thirty_list, key=itemgetter('resolution'))
    thirty_twenty_list = sorted(thirty_twenty_list, key=itemgetter('resolution'))
     
    header = ['group', 'pdb_id', 'uniprot_id', 'chain', 'resolution', 
              'seq_identity', 'seq_identity_percent', 'query_len', 
              'template_len','template_ali_start', 'template_ali_end', 
              'query_sequence', 'alignment', 'template_sequence', 
              'pdb_sequence', 'uniprot_sequence']
     
    with open("alignment_info_v2.csv", 'w') as outfile:
        for dict in ninety_seventy_list:
            if dict["resolution"] > 0:
                for head in header:
                    outfile.write(head)
                    outfile.write(",")
                    if head in dict:
                        outfile.write(str(dict[head]))
                    outfile.write("\n")
                outfile.write("\n")
                outfile.write("\n")
        outfile.write("\n")
        for dict in seventy_fifty_list:
            if dict["resolution"] > 0:
                for head in header:
                    outfile.write(head)
                    outfile.write(",")
                    if head in dict:
                        outfile.write(str(dict[head]))
                    outfile.write("\n")
                outfile.write("\n")
                outfile.write("\n")
        outfile.write("\n")
        for dict in fifty_thirty_list:
            if dict["resolution"] > 0:
                for head in header:
                    outfile.write(head)
                    outfile.write(",")
                    if head in dict:
                        outfile.write(str(dict[head]))
                    outfile.write("\n")
                outfile.write("\n")
                outfile.write("\n")
        outfile.write("\n")
        for dict in thirty_twenty_list:
            if dict["resolution"] > 0:
                for head in header:
                    outfile.write(head)
                    outfile.write(",")
                    if head in dict:
                        outfile.write(str(dict[head]))
                    outfile.write("\n")
                outfile.write("\n")
                outfile.write("\n")
    
    print "Done! :)"
            
            
def get_sequences(file_name):
    seqs = []
    with open(file_name, 'r') as f:
        previous = ""
        query_seq = ""
        template_seq = ""
        template_name = ""
        template_ali_start = []
        template_ali_end = []
        read = False
        count = 0
        for line in f:
            line = re.sub(r"\s+", '\t', line)
            count += 1
            if line.startswith(">"):
                if template_name != "":
                    seqs.append([template_name, 
                                (query_seq, template_seq),
                                (template_ali_start[0], template_ali_end[-1].replace("\n",""))
                                ])
                    print template_ali_start[0] + " " + template_ali_end[-1]
                query_seq = ""
                template_seq = ""
                template_ali_start = []
                template_ali_end = []
                template_name = line.split(";")[0].split("\t")[0].replace(">","")
            if line.startswith("Q\t") and not previous.startswith("Q\t"):
                query = re.sub(r"\s+", '\t', f.next()).split("\t")
                query_seq += query[3]
            if line.startswith("T\t") and not previous.startswith("T\t"):
                template = re.sub(r"\s+", '\t', f.next()).split("\t")
                template_seq += template[3]
                template_ali_start.append(template[2])
                template_ali_end.append(template[4])
            previous = line
        seqs.append([template_name, (query_seq,template_seq)])
    return seqs

def pdb_downloader(query_dict, group):
    if not os.path.exists("pdb_data"):
        os.makedirs("pdb_data")
    group = "pdb_data\\" + group
    if not os.path.exists(group):
        os.makedirs(group)
    query = query_dict.keys()
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = """
    <orgPdbQuery>
    <queryType>org.pdb.query.simple.StructureIdQuery</queryType>
    <description>Simple query for a list of PDB IDs : {0}</description>
    <structureIdList>{0}</structureIdList>
    </orgPdbQuery>
    """.format(query)
    print "query:\n", queryText
    print "querying PDB...\n"
    
    req = urllib2.Request(url, data=queryText)
    f = urllib2.urlopen(req)
    
    result = f.read()
    result = result.split("\n")[:-1]
    print result
    
    pdbl = PDBList()
    PDBlist2=result
    available_pddbs = os.listdir(group)
    
    
    for i in PDBlist2:
        if "pdb" + i.lower() + ".ent" not in available_pddbs:
            pdbl.retrieve_pdb_file(i, file_format="pdb", pdir= group)
        
def pdb_parser(query_dict, group):
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    group = "pdb_data\\" + group
    print "Parsing pdb-files"
    for filename in os.listdir(group):
        structure = parser.get_structure(filename.split(".")[:-1][0][3:], group +"\\"+ filename)
        pdb_record = SeqIO.parse(group +"\\"+ filename, "pdb-seqres")
        pdb_id = structure.id.upper()
        pdb_resolution = structure.header["resolution"]
        try:
            uniprot_id = pdb_record.next().dbxrefs[0]
        except StopIteration:
            uniprot_id = ""
        except IndexError:
            uniprot_id = ""
        try:
            if uniprot_id.split(":")[0] == "UNP":
                uniprot_id = uniprot_id.split(":")[1]
                req = urllib2.Request("https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta")
                f = urllib2.urlopen(req)
                uniprot_sequence = str(SeqIO.parse(StringIO(f.read()), "fasta").next().seq)
            else:
                uniprot_sequence = ""
                uniprot_id = ""
        except urllib2.HTTPError:
            uniprot_sequence = ""
            uniprot_id = ""            
            print("URL", req.get_full_url(), "could not be read.")
        
        pdb_sequences = {}
        
        sequences = ppb.build_peptides(structure)
        
            
        chains = []
        for chain in structure.get_chains():
            chains.append(chain.id)
        chains = set(chains)
        
        if len(sequences) == 0:
            print "Could not extract PDB-sequence for " + pdb_id
            for i, chain in enumerate(chains):
                pdb_sequences[chain] = None
        else:
            for i, chain in enumerate(chains):
                pdb_sequences[chain] = str(sequences[i].get_sequence())
            
        query_dict[pdb_id]["uniprot_id"] = uniprot_id
        query_dict[pdb_id]["resolution"] = pdb_resolution
        query_dict[pdb_id]["pdb_sequence"] = pdb_sequences[query_dict[pdb_id]["chain"]]
        query_dict[pdb_id]["uniprot_sequence"] = uniprot_sequence
        
if __name__ == '__main__':
    main()