import random
from src.dna import Oligonucleotide


class Pool:

    def __init__(self, oligonucleotides: list[Oligonucleotide]):
        self.oligonucleotides = oligonucleotides

    def polymerase_elongation(self):
        i = random.randint(0, len(self.oligonucleotides))

        self.oligonucleotides[i].perform_elongation()

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

    def add_oligonucleotides(self, oligonucleotides):
        self.oligonucleotides += oligonucleotides

    def print(self):
        print("Printing...")
        dnas = ""
        for i, oligo in enumerate(self.oligonucleotides):
            print(str(i) + " of " + str(len(self.oligonucleotides)))

            if oligo.front_left is None:
                dna = ""
                curr = oligo
                while curr is not None:
                    dna += curr.name + ","
                    curr = curr.front_right
                dna += "\n"
                dnas += dna

        print("Writing to file...")
        with open("Output.txt", "w") as text_file:
            text_file.write(dnas)

    def print_strands(self):
        for oligo in self.oligonucleotides:
            first_strand = ""
            second_strand = ""

            if oligo.front_left is None:
                curr = oligo
                while curr is not None:
                    first_strand += curr.name + ":" + "".join(curr.nucleotides) + " "
                    curr = curr.right

                if oligo.front_right is not None:
                    curr = oligo.front_right
                    while curr is not None:
                        second_strand += curr.name + ":" + "".join(curr.nucleotides) + " "
                        curr = curr.right

                if first_strand:
                    print(first_strand)

                if second_strand:
                    print(second_strand)

                print("")
