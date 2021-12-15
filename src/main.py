from tsp import TravelSalesmanProblem


def main():
    problem = TravelSalesmanProblem()
    dna_strands = problem.resolve(
        adj_list={
            '0': ['1', '2'],
            '1': ['2'],
            '2': []
        },
        start='0',
        end='2'
    )


if __name__ == '__main__':
    main()
