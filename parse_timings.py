import argparse
import pstats

valid_sorts = ['calls', 'cumulative', 'cumtime', 'file', 'filename', 'module', 'ncalls', 'pcalls', 'line', 'name', 'nfl', 'stdname', 'time', 'tottime']

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help="file to parse")
parser.add_argument('--sort', type=str, default='cumulative', help="how to sort the results",
					choices=valid_sorts)
parser.add_argument('--count', type=int, default=10, help="total results to display")

args = parser.parse_args()

p = pstats.Stats(args.file)
p.sort_stats(args.sort).print_stats(args.count)
