from typing import List, Dict
from numpy import random

from dna import Oligonucleotide, print_dna_strands
from pool import Pool, Gel


class TravelSalesmanProblem:

    def __init__(self):
        self.multiplication_factor = 0
        self.dna_size = 10
        self.reaction_time = 20
        self.pool = None
        self.n_amplifications = 2
        self.initial_oligonucleotides = None
        self.vertices = None

    def resolve(self, adj_list: Dict[chr, List[chr]], start, end):
        print("Starting...")
        print("Dna Size: " + str(self.dna_size))

        self.vertices, self.initial_oligonucleotides = self.create_oligonucleotides(adj_list, self.dna_size, start, end)

        self.pool = Pool(self.initial_oligonucleotides)
        self.create_random_paths()
        self.amplify_start_and_end_dna(start, end)
        self.get_strands_with_size(len(adj_list.keys()))

    def create_random_paths(self):
        print("Creating random paths...")
        self.pool.annealing(reaction_time=self.reaction_time)
        print_dna_strands(self.pool.get_dna_strands())

    def amplify_start_and_end_dna(self, start, end):
        start_primer = self.create_primer(self.vertices[start])
        end_primer = self.create_primer(Oligonucleotide.complement(self.vertices[end]))

        for i in range(self.n_amplifications):
            print("Denaturing...")
            self.pool.denaturing()
            print_dna_strands(self.pool.get_dna_strands())

            print("Adding primers")
            self.pool.add_oligonucleotides(Oligonucleotide.copy(start_primer) for x in range(2 ^ i))
            self.pool.add_oligonucleotides(Oligonucleotide.copy(end_primer) for x in range(2 ^ i))
            self.pool.annealing_primers(self.reaction_time)
            print_dna_strands(self.pool.get_dna_strands())


            print("Performing polymerase_elongation")

            self.pool.polymerase_chain_reaction()
            print_dna_strands(self.pool.get_dna_strands())

    def get_strands_with_size(self, size):
        print("Placing strands into gel")
        gel = Gel(self.pool.get_dna_strands())

        gel.run()

        print(str(size))
        print_dna_strands(gel.get_strand_with_size(size))

    # Create a list of oligonucleotides representing the edges and the complements
    def create_oligonucleotides(self, adj_list: Dict[chr, List[chr]], dna_size, start, end):
        vertices = {}
        oligonucleotides = []

        # Create Oligonucleotides for each vertex
        for name in adj_list.keys():
            vertices[name] = self.create_random_oligonucleotide(name=name, size=dna_size)

        # Add the DNA edges and complements to the mix
        for vertex_name, paths in adj_list.items():
            complement = Oligonucleotide.complement(vertices[vertex_name])
            complement.name = complement.name + "-co"

            for destination_vertex in paths:
                # Add DNA edge to the mix
                oligonucleotides.append(
                    self.create_edge_oligonucleotide(
                        vertices[vertex_name],
                        vertices[destination_vertex],
                        start,
                        end
                    )
                )

                if destination_vertex == end:
                    end_oligo = Oligonucleotide.complement(vertices[destination_vertex])
                    end_oligo.name += "-co"
                    oligonucleotides.append(end_oligo)

                # For each edge i->j we add i-complement to the mix
                oligonucleotides.append(Oligonucleotide.copy(complement))

        return vertices, self.multiply_oligonucleotides(oligonucleotides)

    def multiply_oligonucleotides(self, oligonucleotides: List[Oligonucleotide]):
        for i in range(self.multiplication_factor):
            oligonucleotides += [Oligonucleotide.copy(item) for item in oligonucleotides]

        return oligonucleotides

    @staticmethod
    def create_primer(oligonucleotide):
        oligonucleotide.is_primer = True
        oligonucleotide.name += "p"
        return oligonucleotide

    @staticmethod
    def create_edge_oligonucleotide(oligo_a, oligo_b, start, end):

        # 3' 10 mers of oligo_a OR 3' 20 mers, if oligo_a == Vin
        left = oligo_a.nucleotides if (oligo_a.name == start) else oligo_a.nucleotides[oligo_a.center_index:]

        # 5' 10 mers of oligo_b OR 5' 20 mers, if oligo_b == Vfin
        right = oligo_b.nucleotides if (oligo_b.name == end) else oligo_b.nucleotides[:oligo_b.center_index]

        return Oligonucleotide(
            nucleotides=left + right,
            center_index=len(left),
            name=oligo_a.name + "->" + oligo_b.name
        )

    @staticmethod
    def create_random_oligonucleotide(name, size):
        return Oligonucleotide(
            nucleotides="".join(random.choice(Oligonucleotide.BASES, size=size)),
            center_index=int(size / 2),
            name=name
        )