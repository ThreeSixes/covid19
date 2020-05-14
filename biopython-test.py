"""
Biopython test script.
"""

import hashlib
import json
import pprint

from Bio import GenBank
from Bio import SeqIO


def extract_feature(feature_record):
    """
    Extract features from GenBank records.
    """

    feature_records = {
        "id": feature_record.id,
        "location": feature_record.location,
        "location_operator": feature_record.location_operator,
        "qualifiers": feature_record.qualifiers,
        "ref": feature_record.ref,
        "ref_db": feature_record.ref_db,
        "strand": feature_record.strand,
        "type": feature_record.type
    }

    return feature_records

def process_genome(file_name, file_format="genbank"):
    """
    Process genome data.
    """

    return_value = {}

    print("Loading genome data and hashing sequences...")

    for index, record in enumerate(SeqIO.parse(file_name, file_format)):
        print(".", end="")

        hasher = hashlib.sha512()
        hasher.update(record.seq.encode("UTF-8"))

        #print("index %i, ID = %s, length %i, with %i features"
        #      % (index, record.id, len(record.seq), len(record.features)))

        record_data = {
            "index"
            "annotations": record.annotations,
            "dbxrefs": record.dbxrefs,
            "description": record.description,
            "features": [],
            "id": record.id,
            "index": index,
            "letter_annotations": record.letter_annotations,
            "lower": record.lower,
            "name": record.name,
            "seq": record.seq,
            "seq_hash": hasher.hexdigest(),
            "upper": record.upper
        }

        for feature in record.features:
            record_data["features"].append(extract_feature(feature))

        return_value.update({record.id: record_data})

    print("")

    return return_value

genomes = process_genome("genome-files/seq_2020-05-13_0236.gb")

genome_hashes = {}

# See which DNA sequnces are indentical.
print("Sorting genomes by DNA sequence hash and writing to file...")
for genome in genomes:
    seq_hash = genomes[genome]['seq_hash']

    if not seq_hash in genome_hashes:
        if 'country' in genomes[genome]['features'][0]['qualifiers']:
            country = genomes[genome]['features'][0]['qualifiers']['country']
        else:
            country = None

        genome_hashes.update({
            genomes[genome]['seq_hash']: {
                "count":   1,
                "ids":     [genomes[genome]['id']],
                "locales": [country]
            }
        })
    else:
        genome_hashes[seq_hash]['count'] += 1
        genome_hashes[seq_hash]['ids'].append(genomes[genome]['id'])
        genome_hashes[seq_hash]['locales'].append(country)

# Write experiment results.
with open("./output/dna_matches.json", "w") as dna_match:
    dna_match.write(json.dumps(genome_hashes))