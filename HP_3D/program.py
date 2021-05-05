from genetic import GeneticAlgorithm
import constant as const
from draw import Draw

class Program:
	def create_new_sequence(self):
		ga = GeneticAlgorithm()
		ga.generation = 0

		Iteration = 1

		ga.initialize()
		ga.evaluate()
		ga.keep_the_best()

		while ga.generation < Iteration:
			ga.generation += 1
			ga.select()
			ga.crossover()
			ga.mutate()
			ga.evaluate()
			ga.elitist()

		ga.report()

		d = Draw()
		d.draw_3D_graphic()
