class Oligonucleotide:
    BASES = ['A', 'T', 'C', 'G']

    def __init__(self, nucleotides: str, center_index: int, name: str):
        self.nucleotides = nucleotides
        self.len = len(nucleotides)
        self.name = name
        self.center_index = center_index
        self.is_primer = False
        self.left = None
        self.right = None
        self.front_left = None
        self.front_right = None
        self.front = None

    def perform_elongation(self):
        if self.front is not None and self.front.is_primer:
            self.front.is_primer = False
            curr = self.right
            while curr is not None:
                curr.anneal_center(Oligonucleotide.complement(curr))
                curr = curr.right

    def denature(self):
        curr = self
        while curr is not None:
            curr.front_left = None
            curr.front_right = None
            curr.front = None
            curr = curr.right

    def anneal(self, other_oligo):
        if self.anneal_left(other_oligo):
            pass
        elif self.anneal_right(other_oligo):
            pass

    def anneal_center(self, other_oligo):
        if self.can_anneal_center(self):
            self.front = other_oligo
            other_oligo.front = self

            other_oligo.left = self.left.front
            other_oligo.right = self.right.front

    def anneal_left(self, other_oligo):
        if self.can_anneal_left(other_oligo):
            # Cross reference them in the diagonal
            self.front_left = other_oligo
            other_oligo.front_right = self

            # Reference the already connected oligos
            other_oligo.right = self.front_right
            if other_oligo.front_left is not None:
                other_oligo.front_left.right = self

    def anneal_right(self, other_oligo):
        if self.can_anneal_right(other_oligo):

            # Cross reference them in the diagonal
            self.front_right = other_oligo
            other_oligo.front_left = self

            # Reference the already connected oligos
            self.right = other_oligo.front_right
            if self.front_left is not None:
                self.front_left.right = other_oligo

    def can_anneal_center(self, other_oligo):
        return self.front_left is not None and self.front_right is not None and self.front is not None and \
            other_oligo.front_left is not None and other_oligo.front_right is not None and other_oligo.front is not None

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
            name=oligonucleotide.name
        )

    @classmethod
    def complement(cls, oligonucleotide):
        return Oligonucleotide(
            nucleotides=cls._inverse(oligonucleotide.nucleotides),
            center_index=oligonucleotide.center_index,
            name=oligonucleotide.name
        )

    @classmethod
    def primer(cls, oligonucleotide):
        oligo = cls.complement(oligonucleotide)
        oligo.is_primer = True
        oligo.name += "-p"
        return oligo
