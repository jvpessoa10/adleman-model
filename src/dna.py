class Oligonucleotide:
    BASES = ['A', 'T', 'C', 'G']

    def __init__(self, nucleotides: str, center_index: int, name: str, is_primer=False):
        self.nucleotides = nucleotides
        self.len = len(nucleotides)
        self.name = name
        self.center_index = center_index
        self.is_primer = is_primer
        self.left = None
        self.right = None
        self.front_left = None
        self.front_right = None
        self.front = None

    def perform_elongation(self):
        if self.is_primer:
            curr = self.front_left if self.front_left is not None \
                else self.front if self.front is not None \
                else self.front_left

            if curr is None:
                return []

            should_go_right = curr.left is None
            curr.denature()  # Removes the Primer
            new_oligos = []
            while curr is not None:
                complement = Oligonucleotide.complement(curr)
                complement.name += "\'"
                curr.anneal_center(complement)
                new_oligos.append(complement)
                curr = curr.right if should_go_right else curr.left


            return new_oligos

    def denature(self):
        self.front_left = None
        self.front_right = None
        self.front = None

    def anneal(self, other_oligo):
        if self.anneal_left(other_oligo):
            pass
        elif self.anneal_right(other_oligo):
            pass
        else:
            self.anneal_center(other_oligo)

    def anneal_center(self, other_oligo):
        if self.can_anneal_center(other_oligo):
            self.front = other_oligo
            other_oligo.front = self

            if self.left is not None:
                other_oligo.left = self.left.front
                if self.left.front is not None:
                    self.left.front.right = other_oligo

            if self.right is not None:
                other_oligo.right = self.right.front

                if self.right.front is not None:
                    self.right.front.left = other_oligo

    def anneal_left(self, other_oligo):
        if self.can_anneal_left(other_oligo):
            # Cross reference them in the diagonal
            self.front_left = other_oligo
            other_oligo.front_right = self

            # Reference the already connected oligos
            other_oligo.right = self.front_right
            if other_oligo.front_left is not None:
                other_oligo.front_left.right = self

            self.left = other_oligo.front_left
            if self.front_right is not None:
                self.front_right.left = other_oligo

    def anneal_right(self, other_oligo):
        if self.can_anneal_right(other_oligo):

            # Cross reference them in the diagonal
            self.front_right = other_oligo
            other_oligo.front_left = self

            # Reference the already connected oligos
            self.right = other_oligo.front_right
            if self.front_left is not None:
                self.front_left.right = other_oligo

            if other_oligo.front_right is not None:
                other_oligo.front_right.left = self

            other_oligo.left = self.front_left

    def can_anneal_center(self, other_oligo):
        return self.front_left is None and self.front_right is None and self.front is None and \
               other_oligo.front_left is None and other_oligo.front_right is None and other_oligo.front is None \
               and self.can_anneal(self.nucleotides, other_oligo.nucleotides)

    def can_anneal_left(self, other_oligo):
        if self.front_left is not None or other_oligo.front_right is not None:
            return False

        return self.can_anneal(self.left_nucleotides(), other_oligo.right_nucleotides()) or \
               (self.can_anneal(self.left_nucleotides(), other_oligo.nucleotides) and other_oligo.front_left is None)

    def can_anneal_right(self, other_oligo):
        if self.front_right is not None or other_oligo.front_left is not None:
            return False

        return self.can_anneal(self.right_nucleotides(), other_oligo.left_nucleotides()) or \
               (self.can_anneal(self.right_nucleotides(), other_oligo.nucleotides) and other_oligo.front_right is None)

    def left_nucleotides(self):
        return self.nucleotides[:self.center_index]

    def right_nucleotides(self):
        return self.nucleotides[self.center_index:]

    @staticmethod
    def can_anneal(strand_a: str, strand_b: str):
        return strand_a == Oligonucleotide._inverse(strand_b)

    @staticmethod
    def _inverse(strand: str) -> str:
        return "".join([Oligonucleotide._inverse_nucleotide(nucleotide) for nucleotide in strand])

    @staticmethod
    def _inverse_nucleotide(nucleotide: chr) -> chr:
        if nucleotide == 'A':
            return 'T'
        elif nucleotide == 'T':
            return 'A'
        elif nucleotide == 'G':
            return 'C'
        elif nucleotide == 'C':
            return 'G'

    @classmethod
    def copy(cls, oligonucleotide):
        return Oligonucleotide(
            nucleotides=oligonucleotide.nucleotides,
            center_index=oligonucleotide.center_index,
            name=oligonucleotide.name,
            is_primer=oligonucleotide.is_primer
        )

    @classmethod
    def complement(cls, oligonucleotide):
        return Oligonucleotide(
            nucleotides=cls._inverse(oligonucleotide.nucleotides),
            center_index=oligonucleotide.center_index,
            name=oligonucleotide.name
        )


def print_dna_strands(strands):
    for (left_strand, right_strand) in strands:
        left_string = ""
        right_string = ""

        curr = left_strand
        while curr is not None:
            left_string += curr.name + ":" + "".join(curr.nucleotides) + " "
            curr = curr.right

        curr = right_strand
        while curr is not None:
            right_string += curr.name + ":" + "".join(curr.nucleotides) + " "
            curr = curr.right

        if left_string:
            print(left_string)
        if right_string:
            print(right_string)
        print("")
