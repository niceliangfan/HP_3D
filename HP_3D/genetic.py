import constant as const
import random
import sqlite3

class Genotype:
	gene = []		# 氨基酸相对于上一个氨基酸的摆放位置
	lower = int(1)		# 氨基酸摆放下限
	upper = int(6)		# 氨基酸摆放得上限
	fitness = int(0)	# 适应度
	rfitness = float(0)	# 相关适应度
	cfitness = float(0)	# 累积适应度

class GeneticAlgorithm:
	const.POPULATIONNUMBER = int(10) # 种群数量(即蛋白质数量)
	const.RESIDUENUMBER = int(360)   # 蛋白质(CB2)的氨基酸数量
	generation = int(0) # 控制遗传的代数
	cur_best = int(0)   # 存放当前所有蛋白质中最优能量值蛋白质的索引
	place = [[[0] * (const.RESIDUENUMBER + 1) for _ in range(const.RESIDUENUMBER + 1) ] for _ in range(const.RESIDUENUMBER + 1)]
	# 三维数组place表示氨基酸在三维空间中的位置
	# place[x, y, z] > 0表示有氨基酸在（x, y, z）的位置上，place[y, x, z] = 0表示没有氨基酸
	# print(place[360, 360, 360], type(place[360, 360, 360]))
	seq = list(range(const.RESIDUENUMBER))
	# 氨基酸相对于上一个氨基酸的位置（六种情况1,2,3,4,5,6）
	protein = list(range(const.RESIDUENUMBER))
	# 用于存放CB2的氨基酸顺序序列
	generation = int(0) # 控制遗传的代数
	cur_best = int(0) # 存放当前所有蛋白质中最优能量值蛋白质的索引
	population = list(range(const.POPULATIONNUMBER)) # 当前种群
	newpopulation = list(range(const.POPULATIONNUMBER)) # 下一代种群
	for i in range(const.POPULATIONNUMBER):
		population[i] = Genotype()
		population[i].gene = list(range(const.RESIDUENUMBER))
		newpopulation[i] = Genotype()
		newpopulation[i].gene = list(range(const.RESIDUENUMBER))

	def createSequence(self, seq, length): # 创建氨基酸的三维空间摆放序列
		result = 0 # 确定氨基酸摆放序列是否创建成功，result=1创建成功，result=0创建失败

		next = 0
		nextp = list(range(10))
		for i in range(10):
			nextp[i] = 0

		for i in range(length + 1):
			for j in range(length + 1):
				for k in range(length + 1):
					self.place[i][j][k] = 0 # 将三维空间中的所有氨基酸清空

		for i in range(length):
			seq[i] = 0 # 将seq(存放氨基酸相对于上一个氨基酸的排放位置)清空

		seq[0] = 9
		# 第一个氨基酸前面没有氨基酸存在，因此其相对于上一个氨基酸的摆放位置没有意义，因此赋值一个与摆放位置无关的数字9
		row = length // 2
		colum = length // 2
		high = length // 2  # 列表的索引必须为整数, type(length/2) = float, type(length//2) = int
		self.place[row][colum][high] = self.protein[0]
		# 第一个氨基酸的位置固定在三维空间的正中央位置

		n = 1
		row = row - 1
		next = 1
		seq[n] = next
		self.place[row][colum][high] = self.protein[n]
		# 第二个氨基酸固定在第一个氨基酸在y轴负方向相差1个单位的位置

		backfrom = 0
		count = 0
		countp = 0
		n = 2
		while n < length:
			next = 0
			countp = self.calculateVacancyNumber(row, colum, high, length, nextp)
			# 在摆放氨基酸之前，检测上一个氨基酸周围的空位数
			if countp > 0: # 有空位置，可以摆放
				next = self.chooseDirection(row, colum, high, length, nextp, countp)
				# 从所有空位中选择出距第一个氨基酸距离最短的空位
				templist = self.dealWithGoNextStep(next, row, colum, high)
				# 根据next修改row、colum、high的值
				row = templist[0]
				colum = templist[1]
				high = templist[2]
				seq[n] = next
				self.place[row][colum][high] = self.protein[n]
			else: # 没有空位，回溯到上一个氨基酸的位置
				self.place[row][colum][high] = 0 # 清空当前位置
				n -= 1
				backfrom = seq[n]
				# 记录氨基酸的当前的摆放位置（1,2,3,4,5,6），下一次摆放时就不再考虑
				templist = self.dealWithGoBackStep(backfrom, row, colum, high)
				# 根据backfrom修改row、colum、high的值
				row = templist[0]
				colum = templist[1]
				high = templist[2]
				count = 1 # 记录回溯次数
				while count <= 40: # 最多进行40次回溯（经检测，回溯40次效果最佳）
					countp = self.calculateVacancyNumber(row, colum, high, length, nextp)
					# 回溯后，再次检测氨基酸附近是否有空位
					nextp[backfrom] = 0
					countp -= 1
					# 刚退回来的位置不予考虑
					if countp > 0: # 有空位置，可以摆放
						next = self.chooseDirection(row, colum, high, length, nextp, countp)
						templist = self.dealWithGoNextStep(next, row, colum, high)
						row = templist[0]
						colum = templist[1]
						high = templist[2]
						seq[n] = next
						self.place[row][colum][high] = self.protein[n]
						break
					else: # 没有空位，继续回溯
						self.place[row][colum][high] = 0
						n -= 1
						backfrom = seq[n]
						templist = self.dealWithGoBackStep(backfrom, row, colum, high)
						row = templist[0]
						colum = templist[1]
						high = templist[2]
						count += 1
				if count > 40: # 回溯40次后仍不成功，则退出序列的创建, 返回值result=0，创建失败
					return result
			n += 1 # 循环(while n < length)+1

		result = 1 # 所有氨基酸均成功摆放，result=1，表示创建成功
		return result

	def checkSequence(self, seq, length): # 检测生成的氨基酸摆放序列是否合法
		for i in range(length + 1):
			for j in range(length + 1):
				for k in range(length + 1):
					self.place[i][j][k] = 0 # 将三维空间中的所有氨基酸清空

		seq[0] = 9
		row = length // 2
		colum = length // 2
		high = length // 2
		self.place[row][colum][high] = self.protein[0]
		# 第一个氨基酸的位置固定在三维空间的正中央位置

		n = 1
		while n < length:
			next = seq[n]
			if next == 1:
				row -= 1
			if next == 2:
				colum += 1
			if next == 3:
				high += 1
			if next == 4:
				row += 1
			if next == 5:
				colum -= 1
			if next == 6:
				high -= 1
			if self.place[row][colum][high] > 0: # 检测到(row, colum, high)位置已存在氨基酸，发生冲突, 氨基酸摆放序列非法, 返回0
				return 0
			self.place[row][colum][high] = self.protein[n]
			n += 1
		return 1

	def initialize(self): # 初始化
		for i in range(const.RESIDUENUMBER + 1):
			for j in range(const.RESIDUENUMBER + 1):
				for k in range(const.RESIDUENUMBER + 1):
					self.place[i][j][k] = 0
		for i in range(const.RESIDUENUMBER):
			self.seq[i] = 0
			self.protein[i] = 0
		for i in range(const.POPULATIONNUMBER):
			for j in range(const.RESIDUENUMBER):
				self.population[i].gene[j] = 0
				self.newpopulation[i].gene[j] = 0
		try:
			f = open('in.txt', 'r')
			string = f.read()
			templist = string.split(' ')
		finally:
			f.close()

		for i in range(const.RESIDUENUMBER):
			self.protein[i] = int(templist[i])

		for i in range(const.POPULATIONNUMBER - 1):
			# print(i)
			while (self.createSequence(self.population[i].gene, const.RESIDUENUMBER) != 1): # 创建氨基酸摆放序列，直至摆放序列合法
				pass

	def randInt(self, low, high): # 生成一个low - high之间的随机数
		# random.seed() # 改变随机数生成器的种子
		result = random.randint(low, high)
		return result

	def randDouble(self): # 生成一个0 - 1之间的随机数
		# random.seed() # 改变随机数生成器的种子
		result = random.random()
		return result

	def evaluate(self): # 评估蛋白质能量值
		for i in range(const.POPULATIONNUMBER - 1):
			resultfit = 0
			if (self.checkSequence(self.population[i].gene, const.RESIDUENUMBER) == 0):
				# 判定氨基酸摆放序列是否创建成功，不成功则重新创建
				while (self.createSequence(self.population[i].gene, const.RESIDUENUMBER) != 1):
					pass
			for j in range(const.RESIDUENUMBER):
				if self.protein[j] != 1: # 判定氨基酸是否为H型氨基酸
					continue
				for k in range(j + 3, const.RESIDUENUMBER):
					if self.protein[k] != 1: # 判定氨基酸是否为H型氨基酸
						continue
					x = 0
					y = 0
					z = 0
					for m in range(j + 1, k + 1):
						if self.population[i].gene[m] == 1:
							y -= 1
						if self.population[i].gene[m] == 2:
							x += 1
						if self.population[i].gene[m] == 3:
							z += 1
						if self.population[i].gene[m] == 4:
							y += 1
						if self.population[i].gene[m] == 5:
							x -= 1
						if self.population[i].gene[m] == 6:
							z -= 1
					if ((x == -1 or x == 1) and y == 0 and z == 0):
						resultfit -= 1
					if ((y == -1 or y == 1) and x == 0 and z == 0):
						resultfit -= 1
					if ((z == -1 or z == 1) and x == 0 and y == 0):
						resultfit -= 1
			self.population[i].fitness = resultfit

	def keep_the_best(self): # 记录当代中能量值最优的蛋白质的摆放序列
		cur_best = 0 # 存放当前最优能量值的蛋白质索引
		for i in range(const.POPULATIONNUMBER - 1):
			if (self.population[i].fitness < self.population[const.POPULATIONNUMBER - 1].fitness):
				cur_best = i
				self.population[const.POPULATIONNUMBER - 1].fitness = self.population[i].fitness
		for i in range(const.RESIDUENUMBER):
			self.population[const.POPULATIONNUMBER - 1].gene[i] = self.population[cur_best].gene[i]

	def elitist(self): # 精英算法，保证一代好于一代
		best = self.population[0].fitness
		worst = self.population[0].fitness
		best_mem = 0
		worst_mem = 0
		for i in range(1, const.POPULATIONNUMBER - 1):
			if (self.population[i].fitness <= best):
				best = self.population[i].fitness
				best_mem = i
			if (self.population[i].fitness >= worst):
				worst = self.population[i].fitness
				worst_mem = i
		# self.population[const.POPULATIONNUMBER - 1].fitness = -150
		# print(self.population[const.POPULATIONNUMBER - 1].fitness)
		# print(best, best_mem)
		# print(worst, worst_mem)
		if (best <= self.population[const.POPULATIONNUMBER - 1].fitness):
			for i in range(const.RESIDUENUMBER):
				self.population[const.POPULATIONNUMBER - 1].gene[i] = self.population[best_mem].gene[i]
			self.population[const.POPULATIONNUMBER - 1].fitness = self.population[best_mem].fitness
		else:
			for i in range(const.RESIDUENUMBER):
				self.population[worst_mem].gene[i] = self.population[const.POPULATIONNUMBER - 1].gene[i]
			self.population[worst_mem].fitness = self.population[const.POPULATIONNUMBER - 1].fitness

	def select(self): # 选择算子
		total = 0
		for i in range(const.POPULATIONNUMBER - 1):
			total += self.population[i].fitness
		# print(total, type(total))

		for i in range(const.POPULATIONNUMBER - 1):
			self.population[i].rfitness = self.population[i].fitness / total
			# print(i, self.population[i].rfitness, type(self.population[i].rfitness))

		self.population[0].cfitness = self.population[0].rfitness
		# print(0, self.population[0].cfitness, type(self.population[0].cfitness))
		for i in range(1, const.POPULATIONNUMBER - 1):
			self.population[i].cfitness = self.population[i - 1].cfitness + self.population[i].rfitness
			# print(i, self.population[i].cfitness, type(self.population[i].cfitness))

		for i in range(const.POPULATIONNUMBER - 1):
			cfitstart = 0
			cfitend = self.population[0].rfitness
			p = self.randDouble()
			for j in range(const.POPULATIONNUMBER - 1):
				if (p >= cfitstart and p < cfitend):
					self.newpopulation[i] = self.population[j + 1]
					for k in range(const.RESIDUENUMBER):
						self.newpopulation[i].gene[k] = self.population[j + 1].gene[k]
					self.newpopulation[i].fitness = self.population[j + 1].fitness
					self.newpopulation[i].rfitness = self.population[j + 1].rfitness
					self.newpopulation[i].cfitness = self.population[j + 1].cfitness
					break
				cfitstart += self.population[j].rfitness
				cfitend = cfitstart + self.population[j + 1].rfitness

		for i in range(const.POPULATIONNUMBER - 1):
			for j in range(const.RESIDUENUMBER):
				self.population[i].gene[j] = self.newpopulation[i].gene[j]
			self.population[i].fitness = self.newpopulation[i].fitness
			self.population[i].rfitness = self.newpopulation[i].rfitness
			self.population[i].cfitness = self.newpopulation[i].cfitness

	def crossover(self): # 交叉算子
		first = 0
		for i in range(const.POPULATIONNUMBER - 1):
			x = self.randDouble()
			if x < 0.1:
				first += 1
				if first % 2 == 0:
					self.xOver(one, i)
				else:
					one = i

	def xOver(self, one, two):
		point = self.randInt(1, const.RESIDUENUMBER)

		for i in range(point, const.RESIDUENUMBER):
			templist = self.swap(self.population[one].gene[i], self.population[two].gene[i])
			self.population[one].gene[i] = templist[2]
			self.population[two].gene[i] = templist[1]
			if (self.checkSequence(self.population[one].gene, const.RESIDUENUMBER) == 1):
				self.population[two].gene[i] = self.population[one].gene[i]
				break
			templist = self.swap(self.population[one].gene[i], self.population[two].gene[i])
			self.population[one].gene[i] = templist[2]
			self.population[two].gene[i] = templist[1]

		for i in range(point, const.RESIDUENUMBER):
			templist = self.swap(self.population[one].gene[i], self.population[two].gene[i])
			self.population[one].gene[i] = templist[2]
			self.population[two].gene[i] = templist[1]
			if (self.checkSequence(self.population[two].gene, const.RESIDUENUMBER) == 1):
				self.population[one].gene[i] = self.population[two].gene[i]
				break
			templist = self.swap(self.population[one].gene[i], self.population[two].gene[i])
			self.population[one].gene[i] = templist[2]
			self.population[two].gene[i] = templist[1]

	def swap(self, x, y):
		result = [0, x, y]
		return result

	def mutate(self): # 变异算子
		result = 0
		for i in range(const.POPULATIONNUMBER - 1):
			for j in range(1, const.RESIDUENUMBER):
				x = self.randDouble()
				if x < 0.000001:
					result = self.existDiagonalVacancy(i, j)
					if result == 1:
						# print('one', 'population', i, 'gene', j)
						break

					result = self.existLTypeVacancy(i,j)
					if result == 1:
						# print('two', 'population', i, 'gene', j)
						break

					result = self.existUTypeVacancy(i,j)
					if result == 1:
						# print('three', 'population', i, 'gene', j)
						break

					result = self.existTrapezoidVacancy(i,j)
					if result == 1:
						# print('four', 'population', i, 'gene', j)
						break

					result = self.randomVariation(i,j)
					if result == 1:
						# print('five', 'population', i, 'gene', j)
						break

	def report(self): # 报告结果
		templist = []
		print('best fitness = %d' %(self.population[const.POPULATIONNUMBER - 1].fitness))
		best_fitness = self.population[const.POPULATIONNUMBER - 1].fitness
		for j in range(const.RESIDUENUMBER):
			print(self.population[const.POPULATIONNUMBER - 1].gene[j], end = ' ')
			templist.append(str(self.population[const.POPULATIONNUMBER - 1].gene[j]))
		print()
		string = ' '.join(templist)

		try:
			f = open('log.txt', 'w') # 打开文件log.txt
			f.write(string) # 将氨基酸摆放序列写入文件log.txt
		finally:
			f.close() # 关闭文件

		conn = sqlite3.connect('data.db')
		# 建立数据库连接
		cursor = conn.cursor()
		# 创建游标cursor
		cursor.execute('insert into result (fitness, seq) values (?, ?)', (best_fitness, string))
		# 将蛋白质中最优能量值及其氨基酸摆放序列插入数据库
		conn.commit()
		# 提交事务
		cursor.close()
		# 关闭游标cursor
		conn.close()
		# 关闭数据库连接

	def calculateVacancyNumber(self, cur_row, cur_colum, cur_high, maxlength, vacancy): # 计算空间中某点周围的空缺个数
		result = 0 # 返回某点周围的空缺个数
		for i in range(len(vacancy)):
			vacancy[i] = 0

		if (cur_row > 0 and self.place[cur_row - 1][cur_colum][cur_high] <= 0):
			result += 1
			vacancy[1] = 1
		if (cur_colum < maxlength and self.place[cur_row][cur_colum + 1][cur_high] <= 0):
			result += 1
			vacancy[2] = 2
		if (cur_high < maxlength and self.place[cur_row][cur_colum][cur_high + 1] <= 0):
			result += 1
			vacancy[3] = 3
		if (cur_row < maxlength and self.place[cur_row + 1][cur_colum][cur_high] <= 0):
			result += 1
			vacancy[4] = 4
		if (cur_colum > 0 and self.place[cur_row][cur_colum -1][cur_high] <= 0):
			result += 1
			vacancy[5] = 5
		if (cur_high > 0 and self.place[cur_row][cur_colum][cur_high - 1] <= 0):
			result += 1
			vacancy[6] = 6
		return result

	def chooseDirection(self, cur_row, cur_colum, cur_high, maxlength, vacancy, nonezeronum):
		# 选择出所有空位中距第一个氨基酸距离最短的空位
		result = 0
		deviationx = cur_row - maxlength // 2
		deviationy = cur_colum - maxlength // 2
		deviationz = cur_high - maxlength // 2
		absmax = self.maxAbsofThree(deviationx, deviationy, deviationz)
		if (absmax == abs(deviationx)):
			if (deviationx >= 0 and vacancy[1] == 1):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 1
				vacancy[8] = 1
			if (deviationx < 0 and vacancy[4] == 4):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 4
				vacancy[8] = 4
		if (absmax == abs(deviationy)):
			if (deviationy >= 0 and vacancy[5] == 5):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 5
				vacancy[8] = 5
			if (deviationy < 0 and vacancy[2] == 2):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 2
				vacancy[8] = 2
		if (absmax == abs(deviationz)):
			if (deviationz >= 0 and vacancy[6] == 6):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 6
				vacancy[8] = 6
			if (deviationz < 0 and vacancy[3] == 3):
				nonezeronum = nonezeronum + 2
				vacancy[7] = 3
				vacancy[8] = 3
		choose = self.randInt(1, nonezeronum)
		j = 0
		for i in range(1, len(vacancy)):
			if vacancy[i] > 0:
				j += 1
			if j == choose:
				result = vacancy[i]
				break
		return result

	def maxAbsofThree(self, x, y, z):
		result = -1
		x = abs(x)
		y = abs(y)
		z = abs(z)
		if (x > y and x > z):
			return x
		if (y > x and y > z):
			return y
		if (z > x and z > y):
			return z
		return result

	def dealWithGoNextStep(self, next, cur_row, cur_colum, cur_high): # 处理氨基酸摆放
		result = [cur_row, cur_colum, cur_high]
		if next == 1:
			result[0] -= 1
		if next == 2:
			result[1] += 1
		if next == 3:
			result[2] += 1
		if next == 4:
			result[0] += 1
		if next == 5:
			result[1] -= 1
		if next == 6:
			result[2] -= 1
		return result

	def dealWithGoBackStep(self, backfrom, cur_row, cur_colum, cur_high): # 处理氨基酸回溯
		result = [cur_row, cur_colum, cur_high]
		if next == 1:
			result[0] += 1
		if next == 2:
			result[1] -= 1
		if next == 3:
			result[2] -= 1
		if next == 4:
			result[0] -= 1
		if next == 5:
			result[1] += 1
		if next == 6:
			result[2] += 1
		return result

	def existDiagonalVacancy(self, pop_num, num): # 变异规则1，存在直角对角空位
		if (num < 1 or num > const.RESIDUENUMBER - 2):
			return 0 # 越界，变异失败
		row = const.RESIDUENUMBER // 2
		colum = const.RESIDUENUMBER // 2
		high = const.RESIDUENUMBER // 2

		for i in range(1, num + 1):
			if (self.population[pop_num].gene[i] == 1):
				row -= 1
			if (self.population[pop_num].gene[i] == 2):
				colum += 1
			if (self.population[pop_num].gene[i] == 3):
				high += 1
			if (self.population[pop_num].gene[i] == 4):
				row += 1
			if (self.population[pop_num].gene[i] == 5):
				colum -= 1
			if (self.population[pop_num].gene[i] == 6):
				high -= 1
		templist = [self.population[pop_num].gene[num], self.population[pop_num].gene[num + 1]]
		# 情况1
		if (templist == [1, 2] and self.place[row + 1][colum + 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum + 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 2
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况2
		if (templist == [2, 4] and self.place[row + 1][colum - 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum - 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 4
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况3
		if (templist == [4, 5] and self.place[row - 1][colum - 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum - 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 5
			self.population[pop_num].gene[num + 1] = 4
			return 1
		# 情况4
		if (templist == [5, 1] and self.place[row - 1][colum + 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum + 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 1
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况5
		if (templist == [3, 2] and self.place[row][colum + 1][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum + 1][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 2
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况6
		if (templist == [2, 6] and self.place[row][colum - 1][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum - 1][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 6
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况7
		if (templist == [6, 5] and self.place[row][colum - 1][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum - 1][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 5
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况8
		if (templist == [5, 3] and self.place[row][colum + 1][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum + 1][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 3
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况9
		if (templist == [3, 1] and self.place[row - 1][colum][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 1
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况10
		if (templist == [1, 6] and self.place[row + 1][colum][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 6
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况11
		if (templist == [6, 4] and self.place[row + 1][colum][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 4
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况12
		if (templist == [4, 3] and self.place[row - 1][colum][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 3
			self.population[pop_num].gene[num + 1] = 4
			return 1
		'''
		前12种情况为num-1、num、num+1三个氨基酸顺时针呈直角摆放(第num个氨基酸为直角顶点), 对角有空位
		情况1-4在xoy平面, 情况5-8在xoz平面, 情况9-12在yoz平面
		'''
		# 情况13
		if (templist == [5, 4] and self.place[row + 1][colum + 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum + 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 4
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况14
		if (templist == [1, 5] and self.place[row + 1][colum - 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum - 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 5
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况15
		if (templist == [2, 1] and self.place[row - 1][colum - 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum - 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 1
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况16
		if (templist == [4, 2] and self.place[row - 1][colum + 1][high] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum + 1][high] = self.protein[num]
			self.population[pop_num].gene[num] = 2
			self.population[pop_num].gene[num + 1] = 4
			return 1
		# 情况17
		if (templist == [5, 6] and self.place[row][colum + 1][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum + 1][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 6
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况18
		if (templist == [3, 5] and self.place[row][colum - 1][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum - 1][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 5
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况19
		if (templist == [2, 3] and self.place[row][colum - 1][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum - 1][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 3
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况20
		if (templist == [6, 2] and self.place[row][colum + 1][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row][colum + 1][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 2
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况21
		if (templist == [4, 6] and self.place[row - 1][colum][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 6
			self.population[pop_num].gene[num + 1] = 4
			return 1
		# 情况22
		if (templist == [3, 4] and self.place[row + 1][colum][high - 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum][high - 1] = self.protein[num]
			self.population[pop_num].gene[num] = 4
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况23
		if (templist == [1, 3] and self.place[row + 1][colum][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row + 1][colum][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 3
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况24
		if (templist == [6, 1] and self.place[row - 1][colum][high + 1] <= 0):
			self.place[row][colum][high] = 0
			self.place[row - 1][colum][high + 1] = self.protein[num]
			self.population[pop_num].gene[num] = 1
			self.population[pop_num].gene[num + 1] = 6
			return 1
		'''
		后12种情况为num-1、num、num+1三个氨基酸逆时针呈直角摆放(第num个氨基酸为直角顶点), 对角有空位
		情况13-16在xoy平面, 情况17-20在xoz平面, 情况21-24在yoz平面
		'''
		return 0 # 上述条件均未满足，变异失败

	def existLTypeVacancy(self, pop_num, num): # 变异规则2，存在L型空位
		if (num < 2 or num > const.RESIDUENUMBER - 2):
			return 0 # 越界，变异失败
		row = const.RESIDUENUMBER // 2
		colum = const.RESIDUENUMBER // 2
		high = const.RESIDUENUMBER // 2

		for i in range(1, num + 1):
			if (self.population[pop_num].gene[i] == 1):
				row -= 1
			if (self.population[pop_num].gene[i] == 2):
				colum += 1
			if (self.population[pop_num].gene[i] == 3):
				high += 1
			if (self.population[pop_num].gene[i] == 4):
				row += 1
			if (self.population[pop_num].gene[i] == 5):
				colum -= 1
			if (self.population[pop_num].gene[i] == 6):
				high -= 1
			if i == num -1:
				row1, colum1, high1 = row, colum, high
			if i == num:
				row2, colum2, high2 = row, colum, high
		templist = [self.population[pop_num].gene[num - 1], self.population[pop_num].gene[num], self.population[pop_num].gene[num + 1]]
		# 情况1
		if (templist == [5, 4, 4] and self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1 + 1][high1] = self.protein[num - 1]
			self.place[row2 + 1][colum2 + 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况2
		if (templist == [1, 5, 5] and self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1 - 1][high1] = self.protein[num - 1]
			self.place[row2 + 1][colum2 - 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况3
		if (templist == [2, 1, 1] and self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1 - 1][high1] = self.protein[num - 1]
			self.place[row2 - 1][colum2 - 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况4
		if (templist == [4, 2, 2] and self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1 + 1][high1] = self.protein[num - 1]
			self.place[row2 - 1][colum2 + 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况5
		if (templist == [5, 6, 6] and self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 + 1][high1 - 1] = self.protein[num - 1]
			self.place[row2][colum2 + 1][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况6
		if (templist == [3, 5, 5] and self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 - 1][high1 - 1] = self.protein[num - 1]
			self.place[row2][colum2 - 1][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况7
		if (templist == [2, 3, 3] and self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 - 1][high1 + 1] = self.protein[num - 1]
			self.place[row2][colum2 - 1][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况8
		if (templist == [6, 2, 2] and self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 + 1][high1 + 1] = self.protein[num - 1]
			self.place[row2][colum2 + 1][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况9
		if (templist == [4, 6, 6] and self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1][high1 - 1] = self.protein[num - 1]
			self.place[row2 - 1][colum2][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 4
			return 1
		# 情况10
		if (templist == [3, 4, 4] and self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1][high1 - 1] = self.protein[num - 1]
			self.place[row2 + 1][colum2][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况11
		if (templist == [1, 3, 3] and self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1][high1 + 1] = self.protein[num - 1]
			self.place[row2 + 1][colum2][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况12
		if (templist == [6, 1, 1] and self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1][high1 + 1] = self.protein[num - 1]
			self.place[row2 - 1][colum2][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 6
			return 1
		'''
		前12种情况为逆时针L型模型的变异方式, 涉及到3个氨基酸：num-1、num、num+1
		情况1-4 xoy平面, 情况5-9 xoz平面, 情况9-12 yoz平面
		'''
		# 情况13
		if (templist == [1, 1, 2] and self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1 + 1][high1] = self.protein[num - 1]
			self.place[row2 + 1][colum2 + 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况14
		if (templist == [2, 2, 4] and self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1 - 1][high1] = self.protein[num - 1]
			self.place[row2 + 1][colum2 - 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况15
		if (templist == [4, 4, 5] and self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1 - 1][high1] = self.protein[num - 1]
			self.place[row2 - 1][colum2 - 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 4
			return 1
		# 情况16
		if (templist == [5, 5, 1] and self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1 + 1][high1] = self.protein[num - 1]
			self.place[row2 - 1][colum2 + 1][high2] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况17
		if (templist == [3, 3, 2] and self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 + 1][high1 - 1] = self.protein[num - 1]
			self.place[row2][colum2 + 1][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况18
		if (templist == [2, 2, 6] and self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 - 1][high1 - 1] = self.protein[num - 1]
			self.place[row2][colum2 - 1][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 2
			return 1
		# 情况19
		if (templist == [6, 6, 5] and self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 - 1][high1 + 1] = self.protein[num - 1]
			self.place[row2][colum2 - 1][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况20
		if (templist == [5, 5, 3] and self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1][colum1 + 1][high1 + 1] = self.protein[num - 1]
			self.place[row2][colum2 + 1][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 5
			return 1
		# 情况21
		if (templist == [3, 3, 1] and self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1][high1 - 1] = self.protein[num - 1]
			self.place[row2 - 1][colum2][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 3
			return 1
		# 情况22
		if (templist == [1, 1, 6] and self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1][high1 - 1] = self.protein[num - 1]
			self.place[row2 + 1][colum2][high2 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 1
			return 1
		# 情况23
		if (templist == [6, 6, 4] and self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 + 1][colum1][high1 + 1] = self.protein[num - 1]
			self.place[row2 + 1][colum2][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 6
			return 1
		# 情况24
		if (templist == [4, 4, 3] and self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
			self.changePartPlace_rule2(row1, colum1, high1, row2, colum2, high2)
			self.place[row1 - 1][colum1][high1 + 1] = self.protein[num - 1]
			self.place[row2 - 1][colum2][high2 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 4
			return 1
		'''
		后12种情况为顺时针L型模型的变异方式, 涉及到3个氨基酸：num-1、num、num+1
		情况13-16 xoy平面, 情况17-20 xoz平面, 情况21-24 yoz平面
		'''
		return 0 # 上述条件均不满足，变异失败

	def changePartPlace_rule2(self, row1, colum1, high1, row2, colum2, high2):
		self.place[row1][colum1][high1] = 0
		self.place[row2][colum2][high2] = 0

	def existUTypeVacancy(self, pop_num, num): # 变异规则3，存在U型空位
		if (num < 2 or num > const.RESIDUENUMBER - 4):
			return 0 # 越界，变异失败
		row = const.RESIDUENUMBER // 2
		colum = const.RESIDUENUMBER // 2
		high = const.RESIDUENUMBER // 2

		for i in range(1, num + 2):
			if (self.population[pop_num].gene[i] == 1):
				row -= 1
			if (self.population[pop_num].gene[i] == 2):
				colum += 1
			if (self.population[pop_num].gene[i] == 3):
				high += 1
			if (self.population[pop_num].gene[i] == 4):
				row += 1
			if (self.population[pop_num].gene[i] == 5):
				colum -= 1
			if (self.population[pop_num].gene[i] == 6):
				high -= 1
			if i == num:
				row1, colum1, high1 = row, colum, high
			if i == num + 1:
				row2, colum2, high2 = row, colum, high

		templist = [self.population[pop_num].gene[num - 1], self.population[pop_num].gene[num], \
		self.population[pop_num].gene[num + 1], self.population[pop_num].gene[num + 2], \
		self.population[pop_num].gene[num + 3]]
		# 情况1
		if (templist == [1, 5, 1, 2, 1]):
			if (self.place[row1][colum1 + 2][high1] <= 0 and self.place[row2][colum2 + 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 2][high1] = self.protein[num]
				self.place[row2][colum2 + 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
			if (self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况2
		if (templist == [2, 1, 2, 4, 2]):
			if (self.place[row1 + 2][colum1][high1] <= 0 and self.place[row2 + 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 2][colum1][high1] = self.protein[num]
				self.place[row2 + 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
			if (self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况3
		if (templist == [4, 2, 4, 5, 4]):
			if (self.place[row1][colum1 - 2][high1] <= 0 and self.place[row2][colum2 - 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 2][high1] = self.protein[num]
				self.place[row2][colum2 - 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况4
		if (templist == [5, 4, 5, 1, 5]):
			if (self.place[row1 - 2][colum1][high1] <= 0 and self.place[row2 - 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 2][colum1][high1] = self.protein[num]
				self.place[row2 - 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况5
		if (templist == [3, 5, 3, 2, 3]):
			if (self.place[row1][colum1 + 2][high1] <= 0 and self.place[row2][colum2 + 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 2][high1] = self.protein[num]
				self.place[row2][colum2 + 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
			if (self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况6
		if (templist == [2, 3, 2, 6, 2]):
			if (self.place[row1][colum1][high1 - 2] <= 0 and self.place[row2][colum2][high2 - 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 - 2] = self.protein[num]
				self.place[row2][colum2][high2 - 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
			if (self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况7
		if (templist == [6, 2, 6, 5, 6]):
			if (self.place[row1][colum1 - 2][high1] <= 0 and self.place[row2][colum2 - 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 2][high1] = self.protein[num]
				self.place[row2][colum2 - 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况8
		if (templist == [5, 6, 5, 3, 5]):
			if (self.place[row1][colum1][high1 + 2] <= 0 and self.place[row2][colum2][high2 + 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 + 2] = self.protein[num]
				self.place[row2][colum2][high2 + 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.popualtion[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.popualtion[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.popualtion[pop_num].gene[num + 2] = 1
				return 1
		# 情况9
		if (templist == [3, 4, 3, 1, 3]):
			if (self.place[row1 - 2][colum1][high1] <= 0 and self.place[row2 - 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 2][colum1][high1] = self.protein[num]
				self.place[row2 - 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况10
		if (templist == [1, 3, 1, 6, 1]):
			if (self.place[row1][colum1][high1 - 2] <= 0 and self.place[row2][colum2][high2 - 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 - 2] = self.protein[num]
				self.place[row2][colum2][high2 - 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
			if (self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况11
		if (templist == [6, 1, 6, 4, 6]):
			if (self.place[row1 + 2][colum1][high1] <= 0 and self.place[row2 + 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 2][colum1][high1] = self.protein[num]
				self.place[row2 + 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
			if (self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况12
		if (templist == [4, 6, 4, 3, 4]):
			if (self.place[row1][colum1][high1 + 2] <= 0 and self.place[row2][colum2][high2 + 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 + 2] = self.protein[num]
				self.place[row2][colum2][high2 + 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		'''
		前12种情况为顺时针U型模型的变异方式, 涉及5个氨基酸：num-1、num、num+1、num+2、num+3
		情况1-4 xoy平面, 情况5-8 xoz平面, 情况9-12 yoz平面
		'''
		# 情况13
		if (templist == [4, 5, 4, 2, 4]):
			if (self.place[row1][colum1 + 2][high1] <= 0 and self.place[row2][colum2 + 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 2][high1] = self.protein[num]
				self.place[row2][colum2 + 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
			if (self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况14
		if (templist == [5, 1, 5, 4, 5]):
			if (self.place[row1 + 2][colum1][high1] <= 0 and self.place[row2 + 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 2][colum1][high1] = self.protein[num]
				self.place[row2 + 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
			if (self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况15
		if (templist == [1, 2, 1, 5, 1]):
			if (self.place[row1][colum1 - 2][high1] <= 0 and self.place[row2][colum2 - 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 2][high1] = self.protein[num]
				self.place[row2][colum2 - 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况16
		if (templist == [2, 4, 2, 1, 2]):
			if (self.place[row1 - 2][colum1][high1] <= 0 and self.place[row2 - 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 2][colum1][high1] = self.protein[num]
				self.place[row2 - 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
		# 情况17
		if (templist == [6, 5, 6, 2, 6]):
			if (self.place[row1][colum1 + 2][high1] <= 0 and self.place[row2][colum2 + 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 2][high1] = self.protein[num]
				self.place[row2][colum2 + 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
			if (self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况18
		if (templist == [5, 3, 5, 6, 5]):
			if (self.place[row1][colum1][high1 - 2] <= 0 and self.place[row2][colum2][high2 - 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 - 2] = self.protein[num]
				self.place[row2][colum2][high2 - 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
			if (self.place[row1 - 1][colum1][high1 - 1] <= 0 and self.place[row2 - 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1][high1 - 1] <= 0 and self.place[row2 + 1][colum2][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 - 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况19
		if (templist == [3, 2, 3, 5, 3]):
			if (self.place[row1][colum1 - 2][high1] <= 0 and self.place[row2][colum2 - 2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 2][high1] = self.protein[num]
				self.place[row2][colum2 - 2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况20
		if (templist == [2, 6, 2, 3, 2]):
			if (self.place[row1][colum1][high1 + 2] <= 0 and self.place[row2][colum2][high2 + 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 + 2] = self.protein[num]
				self.place[row2][colum2][high2 + 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1 - 1][colum1][high1 + 1] <= 0 and self.place[row2 - 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 - 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 + 1][colum1][high1 + 1] <= 0 and self.place[row2 + 1][colum2][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1][high1 + 1] = self.protein[num]
				self.place[row2 + 1][colum2][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
		# 情况21
		if (templist == [6, 4, 6, 1, 6]):
			if (self.place[row1 - 2][colum1][high1] <= 0 and self.place[row2 - 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 2][colum1][high1] = self.protein[num]
				self.place[row2 - 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 1
				self.population[pop_num].gene[num + 2] = 4
				return 1
			if (self.place[row1 - 1][colum1 - 1][high1] <= 0 and self.place[row2 - 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 - 1][colum1 + 1][high1] <= 0 and self.place[row2 - 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 - 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 - 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况22
		if (templist == [4, 3, 4, 6, 4]):
			if (self.place[row1][colum1][high1 - 2] <= 0 and self.place[row2][colum2][high2 - 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 - 2] = self.protein[num]
				self.place[row2][colum2][high2 - 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 6
				self.population[pop_num].gene[num + 2] = 3
				return 1
			if (self.place[row1][colum1 - 1][high1 - 1] <= 0 and self.place[row2][colum2 - 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 + 1][high1 - 1] <= 0 and self.place[row2][colum2 + 1][high2 - 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 - 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 - 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况23
		if (templist == [3, 1, 3, 4, 3]):
			if (self.place[row1 + 2][colum1][high1] <= 0 and self.place[row2 + 2][colum2][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 2][colum1][high1] = self.protein[num]
				self.place[row2 + 2][colum2][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 4
				self.population[pop_num].gene[num + 2] = 1
				return 1
			if (self.place[row1 + 1][colum1 - 1][high1] <= 0 and self.place[row2 + 1][colum2 - 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 - 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 - 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1 + 1][colum1 + 1][high1] <= 0 and self.place[row2 + 1][colum2 + 1][high2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1 + 1][colum1 + 1][high1] = self.protein[num]
				self.place[row2 + 1][colum2 + 1][high2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		# 情况24
		if (templist == [1, 6, 1, 3, 1]):
			if (self.place[row1][colum1][high1 + 2] <= 0 and self.place[row2][colum2][high2 + 2] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1][high1 + 2] = self.protein[num]
				self.place[row2][colum2][high2 + 2] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 3
				self.population[pop_num].gene[num + 2] = 6
				return 1
			if (self.place[row1][colum1 - 1][high1 + 1] <= 0 and self.place[row2][colum2 - 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 - 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 - 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 5
				self.population[pop_num].gene[num + 2] = 2
				return 1
			if (self.place[row1][colum1 + 1][high1 + 1] <= 0 and self.place[row2][colum2 + 1][high2 + 1] <= 0):
				self.changePartPlace_rule3(row1, colum1, high1, row2, colum2, high2)
				self.place[row1][colum1 + 1][high1 + 1] = self.protein[num]
				self.place[row2][colum2 + 1][high2 + 1] = self.protein[num + 1]
				self.population[pop_num].gene[num] = 2
				self.population[pop_num].gene[num + 2] = 5
				return 1
		'''
		后12种情况为逆时针U型模型的变异方式, 涉及5个氨基酸：num-1、num、num+1、num+2、num+3
		情况13-16 xoy平面, 情况17-20 xoz平面, 情况21-24 yoz平面
		'''
		return 0 # 以上条件均不满足，变异失败

	def changePartPlace_rule3(self, row1, colum1, high1, row2, colum2, high2):
		self.place[row1][colum1][high1] = 0
		self.place[row2][colum2][high2] = 0

	def existTrapezoidVacancy(self, pop_num, num): # 变异规则4，存在梯形空位
		if (num < 7 or num > const.RESIDUENUMBER - 2):
			return 0 # 越界，变异失败
		row = const.RESIDUENUMBER // 2
		colum = const.RESIDUENUMBER // 2
		high = const.RESIDUENUMBER // 2

		for i in range(1, num + 1):
			if self.population[pop_num].gene[i] == 1:
				row -= 1
			if self.population[pop_num].gene[i] == 2:
				colum += 1
			if self.population[pop_num].gene[i] == 3:
				high += 1
			if self.population[pop_num].gene[i] == 4:
				row += 1
			if self.population[pop_num].gene[i] == 5:
				colum -= 1
			if self.population[pop_num].gene[i] == 6:
				high -= 1
			if i == num - 3:
				row1, colum1, high1 = row, colum, high
			if i == num - 2:
				row2, colum2, high2 = row, colum, high
			if i == num - 1:
				row3, colum3, high3 = row, colum, high
			if i == num:
				row4, colum4, high4 = row, colum, high

		templist = [self.population[pop_num].gene[num - 6], self.population[pop_num].gene[num - 5], \
		self.population[pop_num].gene[num - 4], self.population[pop_num].gene[num - 3], \
		self.population[pop_num].gene[num - 2], self.population[pop_num].gene[num - 1], \
		self.population[pop_num].gene[num], self.population[pop_num].gene[num + 1]]
		# 情况1
		if (templist == [1, 1, 2, 1, 2, 4, 4, 4] and \
		self.place[row3 + 1][colum3 + 1][high3] <= 0 and self.place[row4 + 1][colum4 + 1][high4] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 + 1][colum3 + 1][high3] = self.protein[num - 1]
			self.place[row4 + 1][colum4 + 1][high4] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 2
			self.population[pop_num].gene[num - 2] = 4
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 5
			return 1

		# 情况2
		if (templist == [2, 2, 4, 2, 4, 5, 5, 5] and \
		self.place[row3 + 1][colum3 - 1][high3] <= 0 and self.place[row4 + 1][colum4 - 1][high4] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 + 1][colum3 - 1][high3] = self.protein[num - 1]
			self.place[row4 + 1][colum4 - 1][high4] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 4
			self.population[pop_num].gene[num - 2] = 5
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 1
			return 1

		# 情况3
		if (templist == [4, 4, 5, 4, 5, 1, 1, 1] and \
		self.place[row3 - 1][colum3 - 1][high3] <= 0 and self.place[row4 - 1][colum4 - 1][high4] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 - 1][colum3 - 1][high3] = self.protein[num - 1]
			self.place[row4 - 1][colum4 - 1][high4] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 5
			self.population[pop_num].gene[num - 2] = 1
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 2
			return 1

		# 情况4
		if (templist == [5, 5, 1, 5, 1, 2, 2, 2] and \
		self.place[row3 - 1][colum3 + 1][high3] <= 0 and self.place[row4 - 1][colum4 + 1][high4] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 - 1][colum3 + 1][high3] = self.protein[num - 1]
			self.place[row4 - 1][colum4 + 1][high4] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 1
			self.population[pop_num].gene[num - 2] = 2
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 4
			return 1

		# 情况5
		if (templist == [3, 3, 2, 3, 2, 6, 6, 6] and \
		self.place[row3][colum3 + 1][high3 - 1] <= 0 and self.place[row4][colum4 + 1][high4 - 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3][colum3 + 1][high3 - 1] = self.protein[num - 1]
			self.place[row4][colum4 + 1][high4 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 2
			self.population[pop_num].gene[num - 2] = 6
			self.population[pop_num].gene[num - 1] = 2
			self.population[pop_num].gene[num + 1] = 5
			return 1

		# 情况6
		if (templist == [2, 2, 6, 2, 6, 5, 5, 5] and \
		self.place[row3][colum3 - 1][high3 - 1] <= 0 and self.place[row4][colum4 - 1][high4 - 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3][colum3 - 1][high3 - 1] = self.protein[num - 1]
			self.place[row4][colum4 - 1][high4 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 6
			self.population[pop_num].gene[num - 2] = 5
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 3
			return 1

		# 情况7
		if (templist == [6, 6, 5, 6, 5, 3, 3, 3] and \
		self.place[row3][colum3 - 1][high3 + 1] <= 0 and self.place[row4][colum4 - 1][high4 + 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3][colum3 - 1][high3 + 1] = self.protein[num - 1]
			self.place[row4][colum4 - 1][high4 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 5
			self.population[pop_num].gene[num - 2] = 3
			self.population[pop_num].gene[num - 1] = 5
			self.population[pop_num].gene[num + 1] = 2
			return 1

		# 情况8
		if (templist == [5, 5, 3, 5, 3, 2, 2, 2] and \
		self.place[row3][colum3 + 1][high3 + 1] <= 0 and self.place[row4][colum4 + 1][high4 + 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3][colum3 + 1][high3 + 1] = self.protein[num - 1]
			self.place[row4][colum4 + 1][high4 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 3
			self.population[pop_num].gene[num - 2] = 2
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 6
			return 1

		# 情况9
		if (templist == [3, 3, 1, 3, 1, 6, 6, 6] and \
		self.place[row3 - 1][colum3][high3 - 1] <= 0 and self.place[row4 - 1][colum4][high4 - 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 - 1][colum3][high3 - 1] = self.protein[num - 1]
			self.place[row4 - 1][colum4][high4 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 1
			self.population[pop_num].gene[num - 2] = 6
			self.population[pop_num].gene[num - 1] = 1
			self.population[pop_num].gene[num + 1] = 4
			return 1

		# 情况10
		if (templist == [1, 1, 6, 1, 6, 4, 4, 4] and \
		self.place[row3 + 1][colum3][high3 - 1] <= 0 and self.place[row4 + 1][colum4][high4 - 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 + 1][colum3][high3 - 1] = self.protein[num - 1]
			self.place[row4 + 1][colum4][high4 - 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 6
			self.population[pop_num].gene[num - 2] = 4
			self.population[pop_num].gene[num - 1] = 6
			self.population[pop_num].gene[num + 1] = 3
			return 1

		# 情况11
		if (templist == [6, 6, 4, 6, 4, 3, 3, 3] and \
		self.place[row3 + 1][colum3][high3 + 1] <= 0 and self.place[row4 + 1][colum4][high4 + 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 + 1][colum3][high3 + 1] = self.protein[num - 1]
			self.place[row4 + 1][colum4][high4 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 4
			self.population[pop_num].gene[num - 2] = 3
			self.population[pop_num].gene[num - 1] = 4
			self.population[pop_num].gene[num + 1] = 1
			return 1

		# 情况12
		if (templist == [4, 4, 3, 4, 3, 1, 1, 1] and \
		self.place[row3 - 1][colum3][high3 + 1] <= 0 and self.place[row4 - 1][colum4][high4 + 1] <= 0):
			self.changePartPlace_rule4(row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4)
			self.place[row3 - 1][colum3][high3 + 1] = self.protein[num - 1]
			self.place[row4 - 1][colum4][high4 + 1] = self.protein[num]
			self.population[pop_num].gene[num - 3] = 3
			self.population[pop_num].gene[num - 2] = 1
			self.population[pop_num].gene[num - 1] = 3
			self.population[pop_num].gene[num + 1] = 6
			return 1
		'''
		梯形空位不存在顺时针、逆时针问题
		情况1-4 xoy平面, 情况5-8 xoz平面, 情况9-12 yoz平面
		'''
		return 0 # 以上条件均不满足，变异失败

	def changePartPlace_rule4(self, row1, colum1, high1, row2, colum2, high2, row3, colum3, high3, row4, colum4, high4):
		self.place[row1][colum1][high1] = 0
		self.place[row2][colum2][high2] = 0
		self.place[row3][colum3][high3] = 0
		self.place[row4][colum4][high4] = 0
		self.place[row3][colum3][high3] = self.protein[num - 3]
		self.place[row4][colum4][high4] = self.protein[num - 2]

	def randomVariation(self, pop_num, num): # 变异规则5，氨基酸摆放位置随机变异
		if (num < 0 or num > const.RESIDUENUMBER - 1):
			return 0 # 越界，变异失败
		for i in range(num, const.RESIDUENUMBER):
			lbound = self.population[pop_num].lower
			hbound = self.population[pop_num].upper
			old_gene = self.population[pop_num].gene[i]
			self.population[pop_num].gene[i] = self.population[pop_num].gene[i-1] + 4 + self.randInt(lbound, hbound)
			self.population[pop_num].gene[i] = self.population[pop_num].gene[i] % 6
			if self.population[pop_num].gene[i] == 0:
				self.population[pop_num].gene[i] = 6
			if (self.checkSequence(self.population[pop_num].gene, const.RESIDUENUMBER) == 1):
				# print('five', 'population', pop_num, 'gene', i)
				return 1
			self.population[pop_num].gene[i] = old_gene
		return 0 # 变异失败
