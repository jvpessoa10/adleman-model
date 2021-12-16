from src.dna import Oligonucleotide


class Pool:

    def __init__(self, oligonucleotides: list[Oligonucleotide]):
        self.oligonucleotides = oligonucleotides

    def denaturing(self):
        # i = random.randint(0, len(self.oligonucleotides))
        # self.oligonucleotides[i].denature()
        for oligo1 in self.oligonucleotides:
            oligo1.denature()

    def annealing(self, reaction_time):
        # start = time.time()
        # while time.time() - start < reaction_time:
        #     i = random.randint(0, len(self.oligonucleotides))
        #     j = random.randint(0, len(self.oligonucleotides))
        #
        #     self.oligonucleotides[i].anneal(self.oligonucleotides[j])
        for oligo1 in self.oligonucleotides:
            for oligo2 in self.oligonucleotides:
                if oligo1 != oligo2:
                    oligo1.anneal(oligo2)

    def annealing_primers(self, reaction_time):
        for oligo in self.oligonucleotides:
            for oligo2 in self.oligonucleotides:
                if oligo != oligo2 and oligo2.is_primer:
                    oligo.anneal(oligo2)

        print("Done")

    def add_oligonucleotides(self, oligonucleotides):
        self.oligonucleotides += oligonucleotides

    def polymerase_chain_reaction(self):
        new_oligos = []
        for oligo in self.oligonucleotides:
            if oligo.is_primer:
                new_oligos += oligo.perform_elongation()

        self.oligonucleotides += new_oligos

    def _get_strands(self):
        return [oligo for oligo in self.oligonucleotides if oligo.left is None]

    def get_dna_strands(self):
        visited = []
        strands = []
        for oligo in self._get_strands():

            if oligo in visited or oligo.front_left in visited or oligo.front in visited:
                continue

            visited.append(oligo)
            curr = oligo
            dna = ()
            while curr is not None:
                dna = dna + (curr,)
                if curr.front_left is not None:
                    strands.append((curr, curr.front_left))
                    visited.append(curr.front_left)
                elif curr.front is not None:
                    strands.append((curr, curr.front))
                    visited.append(curr.front)
                elif curr.front_right is not None:
                    strands.append((curr, curr.front_right))
                    visited.append(curr.front_right)

                curr = curr.right
        return strands


class Gel:
    def __init__(self, strands: list[(Oligonucleotide, Oligonucleotide)]):
        self.strands = strands
        self.strands_with_distance = []

    def run(self):
        for strand_left, strand_right in self.strands:
            size_left = self._get_strand_size(strand_left)
            size_right = self._get_strand_size(strand_right)

            size = size_left if size_left > size_right else size_right

            self.strands_with_distance.append((strand_left, strand_right, size))

    def get_strand_with_size(self, size):
        return [(left, right) for (left, right, distance) in self.strands_with_distance if distance == size]

    def _get_strand_size(self, strand):
        size = 0
        curr = strand
        while curr is not None:
            size += 1
            curr = curr.right

        return size
