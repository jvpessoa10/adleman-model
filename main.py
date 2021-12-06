from dna import TSP, Oligonucleotide


def test_dna():
    first = Oligonucleotide.from_list(["A", "A", "G", "G"], name="First")
    second = Oligonucleotide.from_list(["T", "T"], name="Second")

    third = Oligonucleotide.from_list(["C", "C"], name="Third")

    first.annealing(second)
    first.annealing(third)

    print(first)


def main():
    problem = TSP()
    problem.resolve(
        {
            'A': ['B'],
            'B': ['C'],
            'C': []
        },
        'A',
        'B'
    )

if __name__ == '__main__':
    main()
