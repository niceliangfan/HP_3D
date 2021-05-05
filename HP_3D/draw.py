import tkinter
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sqlite3

class Draw:
	def draw_3D_graphic(self):
		# 读取文件in.txt中蛋白质的氨基酸顺序序列
		try:
			f = open('in.txt', 'r')
			string = f.read()
			# 通过空格分离字符串
			amino_acidSequence = string.split(' ')
		finally:
			f.close()
		# 读取文件log.txt中蛋白质的氨基酸摆放序列
		try:
			f = open('log.txt', 'r')
			string = f.read()
			# 通过空格分离字符串
			amino_acidPlacement = string.split(' ')
		finally:
			f.close()

		fig = plt.figure()
		ax = Axes3D(fig)

		x = 180
		y = 180
		z = 180
		ax.scatter(x, y, z, c = 'r')
		# 绘制第一个氨基酸, 固定在三维空间(360*360*360)中的正中央

		for i in range(1, 360):
			x1, y1, z1 = x, y, z
			# H型氨基酸(1)为红色, P型氨基酸(2)为绿色
			if amino_acidSequence[i] == '1':
				color = 'r'
			else:
				color = 'g'
			# 氨基酸摆放规则, 1:x--; 2:y++; 3:z++; 4:x++; 5:y--; 6:z--
			if amino_acidPlacement[i] == '1':
				x -= 1
			if amino_acidPlacement[i] == '2':
				y += 1
			if amino_acidPlacement[i] == '3':
				z += 1
			if amino_acidPlacement[i] == '4':
				x += 1
			if amino_acidPlacement[i] == '5':
				y -= 1
			if amino_acidPlacement[i] == '6':
				z -= 1
			x2, y2, z2 = x, y, z
			x0 = [x1, x2]
			y0 = [y1, y2]
			z0 = [z1, z2]
			ax.scatter(x, y, z, c = color) # 绘制氨基酸
			ax.plot(x0, y0, z0, c='c') # 连接两个氨基酸

		# 绘制XYZ坐标轴
		ax.set_xlabel('X', fontdict={'size': 15, 'color': 'black'})
		ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'black'})
		ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'black'})
		# 显示结果
		plt.show()

	def show_best_sequence(self):
		# 读取文件in.txt中蛋白质的的氨基酸顺序序列
		try:
			f = open('in.txt', 'r')
			string = f.read()
			# 通过空格分离字符串
			amino_acidSequence = string.split(' ')
		finally:
			f.close()

		try:
			# 从数据库中读取蛋白质的最优摆放序列
			conn = sqlite3.connect('data.db')
			cursor = conn.cursor()
			cursor.execute('select  * from result_best_fitness_view')
			values = cursor.fetchall()
			string = values[0][1]
			# 通过空格分离字符串
			amino_acidPlacement = string.split(' ')
		finally:
			cursor.close()
			conn.close()

		fig = plt.figure()
		ax = Axes3D(fig)

		x = 180
		y = 180
		z = 180
		ax.scatter(x, y, z, c = 'r')
		# 绘制第一个氨基酸, 固定在三维空间(360*360*360)中的正中央

		for i in range(1, 360):
			x1, y1, z1 = x, y, z
			# H型氨基酸(1)为红色, P型氨基酸(2)为绿色
			if amino_acidSequence[i] == '1':
				color = 'r'
			else:
				color = 'g'
			# 氨基酸摆放规则, 1:x--; 2:y++; 3:z++; 4:x++; 5:y--; 6:z--
			if amino_acidPlacement[i] == '1':
				x -= 1
			if amino_acidPlacement[i] == '2':
				y += 1
			if amino_acidPlacement[i] == '3':
				z += 1
			if amino_acidPlacement[i] == '4':
				x += 1
			if amino_acidPlacement[i] == '5':
				y -= 1
			if amino_acidPlacement[i] == '6':
				z -= 1
			x2, y2, z2 = x, y, z
			x0 = [x1, x2]
			y0 = [y1, y2]
			z0 = [z1, z2]
			ax.scatter(x, y, z, c = color) # 绘制氨基酸
			ax.plot(x0, y0, z0, c='c') # 连接两个氨基酸

		# 绘制XYZ坐标轴
		ax.set_xlabel('X', fontdict={'size': 15, 'color': 'black'})
		ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'black'})
		ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'black'})
		# 显示结果
		plt.show()
