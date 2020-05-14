"""
Biopython test script.
"""

import hashlib
import json
import pprint

from Bio import SeqIO


def process_genome(file_name, file_type="genbank"):
    """
    """

    return_value = {}

    record_dict = SeqIO.to_dict(SeqIO.parse(file_name, file_type))

    print("Found %s records in %s." %(len(record_dict), file_name))
    for record in record_dict:
        dna_hash = hashlib.sha256()
        dna_hash.update(record_dict[record].seq.encode('UTF-8'))

        record_content = {
            "annotations": record_dict[record].annotations,
            "dna_hash":    dna_hash.hexdigest(),
            "id":          record_dict[record].id,
            "description": record_dict[record].description,
            "name":        record_dict[record].name,
            "seq":         record_dict[record].seq
        }

        return_value.update({record: record_content})

    return return_value

genomes = process_genome("genome-files/seq_2020-05-13_0236.gb")
genome_hashes = {}

for genome in genomes:
    dna_hash = genomes[genome]['dna_hash']

    if not dna_hash in genome_hashes:
        genome_hashes.update({
            genomes[genome]['dna_hash']: {
                "count": 1,
                "ids": [genomes[genome]['id']]
            }
        })
    else:
        genome_hashes[dna_hash]['count'] += 1
        genome_hashes[dna_hash]['ids'].append(genomes[genome]['id'])

pprint.pprint(genome_hashes)