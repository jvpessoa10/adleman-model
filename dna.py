from typing import List, Dict
from numpy import random

BASES = ['A', 'T', 'G', 'C']


class Nucleotide:

    def __init__(self, base: chr):
        self.base = base
        self.adjacent_nucleotide = None
        self.oligonucleotide = None

    def bind(self, nucleotide):
        if self.can_bind(nucleotide):
            self.adjacent_nucleotide = nucleotide
            nucleotide.adjacent_nucleotide = self
            return True

        return False

    def can_bind(self, nucleotide):
        return self.can_bind_base(nucleotide.base) or \
               (self.adjacent_nucleotide is None and nucleotide.adjacent_nucleotide is None)

    def can_bind_base(self, base):
        return self.inverse() == base

    def inverse(self):
        if self.base == 'A':
            return 'T'
        elif self.base == 'T':
            return 'A'
        elif self.base == 'G':
            return 'C'
        elif self.base == 'C':
            return 'G'


class Oligonucleotide:

    def __init__(self, nucleotides: List[Nucleotide], name: str):
        self.nucleotides = nucleotides
        self.name = name

        for nucleotide in nucleotides:
            nucleotide.oligonucleotide = self

    def annealing(self, oligonucleotide):
        valid_position = self._check_valid_space(oligonucleotide)

        if valid_position != -1:
            self._perform_annealing(oligonucleotide, valid_position)

    def _perform_annealing(self, oligonucleotide, pos):
        for i, nucleotide in enumerate(self.nucleotides[pos:]):
            if i < len(oligonucleotide.nucleotides):
                nucleotide.bind(oligonucleotide.nucleotides[i])

    def _check_valid_space(self, oligonucleotide, pos=0):
        if pos == len(self.nucleotides) - 1:
            return -1

        strand = oligonucleotide.nucleotides

        for i, nucleotide in enumerate(self.nucleotides[pos:]):
            if i >= len(strand):
                return pos

            if not nucleotide.can_bind(strand[i]):
                return self._check_valid_space(oligonucleotide, pos + 1)

        return pos

    def complement(self):
        return Oligonucleotide.from_list(
            [nucleotide.inverse() for nucleotide in self.nucleotides],
            name=self.name + "-complement"
        )

    @classmethod
    def from_list(cls, bases: list[chr], name: str):
        return Oligonucleotide(
            [Nucleotide(base) for base in bases],
            name=name
        )

    @classmethod
    def random(cls, size, name: str):
        return Oligonucleotide.from_list(
            random.choice(BASES, size=size).tolist(),
            name=name
        )


class DNAGraph:
    def __init__(self):
        self.vertices: List[Oligonucleotide] = []
        self.paths: List[Oligonucleotide] = []

    def create(self, adj_list: Dict[chr, List[chr]], dna_size, start, end):
        vertices = {}
        paths = []

        # Create DNA for each vertex
        for vertex in adj_list.keys():
            vertices[vertex] = self._create_random_dna_vertex(dna_size, name=vertex)

        for node, adj_nodes in adj_list.items():
            for adj_node in adj_nodes:
                path = self._create_dna_path(
                    vertices[node],
                    vertices[adj_node],
                    dna_size=dna_size,
                    start=node == start,
                    end=adj_node == end
                )
                paths.append(path)

        self.vertices = list(vertices.values())
        self.paths = paths

    def _create_dna_path(self, a: Oligonucleotide, b: Oligonucleotide, dna_size, start, end):
        if start:
            left_part = a.nucleotides
        else:
            left_part = a.nucleotides[int(dna_size / 2):]

        if end:
            right_part = b.nucleotides
        else:
            right_part = b.nucleotides[:int(dna_size / 2)]

        return self._create_dna_vertex(
            left_part + right_part,
            a.name + "->" + b.name
        )

    @staticmethod
    def _create_dna_vertex(bases: list[Nucleotide], name):
        return Oligonucleotide(bases, name)

    @staticmethod
    def _create_random_dna_vertex(dna_size, name):
        return Oligonucleotide.random(dna_size, name=name)


class Pool:
    def __init__(self, oligonucleotides: list[Oligonucleotide]):
        self.oligonucleotides = oligonucleotides

    def start_reaction(self):
        for i in self.oligonucleotides:
            for j in self.oligonucleotides:
                if i != j:
                    i.annealing(j)

    def print(self):
        for oligonucleotide in self.oligonucleotides:
            print("Name:" + oligonucleotide.name)
            print("Nucleotides:")
            for nucleotide in oligonucleotide.nucleotides:
                msg = oligonucleotide.name + " " + nucleotide.base + "-"
                if nucleotide.adjacent_nucleotide is not None:
                    msg += nucleotide.adjacent_nucleotide.base
                    msg += " " + nucleotide.adjacent_nucleotide.oligonucleotide.name
                print(msg)
            print("")


class TSP:
    DNA_SIZE = 4

    def __init__(self):
        self.graph = None
        self.pool = None

    def resolve(self, adj_list: Dict[chr, List[chr]], start, end):
        self.graph = DNAGraph()
        self.graph.create(adj_list, TSP.DNA_SIZE, start, end)

        self.pool = Pool(self.graph.paths + [item.complement() for item in self.graph.vertices])
        self.pool.start_reaction()
        self.pool.print()
