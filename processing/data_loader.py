import random
import os

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez

from abc import ABC, abstractmethod


class DataLoader(ABC):
    @abstractmethod
    def load(self, *args):
        pass


class RandomSequenceLoader(DataLoader):
    def __init__(self):
        self.alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def load(self, length):
        return Seq(''.join(random.choices(self.alphabet, k=length)))
    

class ApiSequenceLoader(DataLoader):
    def __init__(self):
        self.email = 'example@gmail.com'
        Entrez.email = self.email

    def load(self, id_list, type):
        ids = ",".join(id_list)
        handle = Entrez.efetch(db=type, id=ids, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))

        handle.close()

        self.get_sequences_by_id(id_list, type)
        return sequences
    
    def get_sequences_by_id(self, ids, db, folder="sequences"):
        routes = []
        for id in ids:
            print(f"Downloading sequence {id} from {db} database")
            with Entrez.efetch(db=db, id=id, rettype="fasta", retmode="text") as fetch:
                content = fetch.read()
            
            first_line = content.split("\n")[0]
            unique_id = first_line.split()[0][1:]
            filename = unique_id.replace(".", "_") + ".fasta"
            file_route = os.path.join(folder, filename)
            
            with open(file_route, "w") as file:
                file.write(content)
                print(f"File {filename} created")

            routes.append(file_route)
        
        return routes


class FastaDataLoader(DataLoader):
    def __init__(self):
        pass

    def load(self, filename):
        lines = []
        with open(filename, "r") as file:
            lines = file.readlines()

        return ''.join(lines[1:]).replace("\n", "")
    

class DataLoaderFactory:
    @staticmethod
    def get_loader(loader_type):
        if loader_type == "random":
            return RandomSequenceLoader()
        
        elif loader_type == "api":
            return ApiSequenceLoader()
        
        elif loader_type == "fasta":
            return FastaDataLoader()