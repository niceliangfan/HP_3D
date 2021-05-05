from tkinter import *
from sqlite3 import *
from tkinter.messagebox import *
from draw import *
from program import *

class OperateFrame(Frame): # 继承Frame类
	def __init__(self,master=None):
		Frame.__init__(self,master)
		self.root = master #定义内部变量root
		self.createPage()

	def createPage(self):
		Label(self, text = '操作界面', font = ('黑体', 15), fg = 'blue').grid(row = 0, column = 1)
		Button(self, text = '创建新序列', command = self.createNewSequence).grid(row = 2, column = 0)
		Button(self, text = '显示最优序列', command = self.showBestSequence).grid(row = 2, column = 2)

	def createNewSequence(self):
		answer = askokcancel(title = '确认', message = '生成新序列的时间较长,是否确认生成?')
		if answer:
			p = Program()
			p.create_new_sequence()

	def showBestSequence(self):
		d = Draw()
		d.show_best_sequence()

class UserFrame(Frame): # 继承Frame类
	def __init__(self, master):
		Frame.__init__(self, master)
		self.root = master
		self.createPage()
	def createPage(self):
		
		self.listbox = Listbox(self, height = 12)
		self.listbox.pack(side = LEFT)

		Label(self, text = '账号:').pack(fill = X)
		self.username_entry = Entry(self)
		self.username_entry.pack()

		Label(self, text = '密码:').pack(fill = X)
		self.password_entry = Entry(self)
		self.password_entry.pack()

		show_button = Button(self, text = '显示所有用户', command = self.show)
		show_button.pack(fill = X)

		add_button = Button(self, text = '增加用户信息', command = self.add)
		add_button.pack(fill = X)

		update_button = Button(self, text = '更新用户信息', command = self.update)
		update_button.pack(fill = X)

		delete_button = Button(self, text = '删除用户信息', command = self.delete)
		delete_button.pack(fill = X)

	def show(self):
		self.listbox.delete(0, END)
		conn = connect('data.db')
		cursor = conn.cursor()
		cursor.execute('select uid, uname, upassword from user')
		values = cursor.fetchall()
		cursor.close()
		conn.close()
		for value in values:
			string = str(value[1]) + ' ' + value[2]
			self.listbox.insert(END, string)

	def add(self):
		if self.username_entry.get() != '' and self.password_entry.get() != '':
			self.listbox.insert(self.listbox.size(), self.username_entry.get() + ' ' + self.password_entry.get())
			conn = connect('data.db')
			cursor = conn.cursor()
			cursor.execute('insert into user (uname, upassword) values(?, ?)', (self.username_entry.get(), self.password_entry.get()))
			conn.commit()
			cursor.close()
			conn.close()
		else:
			showinfo(title='提示',message='账号或密码不能为空！')

	def update(self):
		if self.listbox.curselection() == ():
			showinfo(title='提示',message='未选中要修改信息的用户！')
		else:	
			if (self.username_entry.get() != '' or self.password_entry.get() != ''):
				selected = self.listbox.curselection()[0]
				conn = connect('data.db')
				cursor = conn.cursor()
				cursor.execute('select * from user')
				values = cursor.fetchall()
		
				self.listbox.delete(selected)
				if self.username_entry.get() != '' and self.password_entry.get() != '':
					self.listbox.insert(selected, self.username_entry.get() + ' ' + self.password_entry.get())
					cursor.execute('update user set uname = ?, upassword = ? where uid = ?', (self.username_entry.get(), self.password_entry.get(), int(values[int(selected)][0])))
				if self.username_entry.get() != '' and self.password_entry.get() == '':
					self.listbox.insert(selected, self.username_entry.get() + ' ' + str(values[int(selected)][2]))
					cursor.execute('update user set uname = ? where uid = ?', (self.username_entry.get(), int(values[int(selected)][0])))
				if self.username_entry.get() == '' and self.password_entry.get() != '':
					self.listbox.insert(selected, str(values[int(selected)][1]) + ' ' + self.password_entry.get())
					cursor.execute('update user set upassword = ? where uid = ?', (self.password_entry.get(), int(values[int(selected)][0])))

				conn.commit()
				cursor.close()
				conn.close()
			else:
				showinfo(title='提示',message='请输入要修改的信息！')

	def delete(self):
		if self.listbox.curselection() != ():
			selected = self.listbox.curselection()[0]
			self.listbox.delete(selected)
			conn = connect('data.db')
			cursor = conn.cursor()
			cursor.execute('select * from user')
			values = cursor.fetchall()
			cursor.execute('delete from user where uid = ?', (int(values[int(selected)][0]),))
			conn.commit()
			cursor.close()
			conn.close()
		else:
			showinfo(title = '提示', message = '请选择要删除的用户')

class SoftWareInfoFrame(Frame): # 继承Frame类
	def __init__(self,master):
		Frame.__init__(self,master)
		self.root = master #定义内部变量root
		self.createPage()

	def createPage(self):
		Label(self,text='该系统使用python语言编写\n前端系统使用tkinter和matplotlib模块\n后台系统基于HP格点模型，使用遗传算法进行优化', justify = LEFT, font = ('宋体', 15)).pack()

class DeveloperInfoFrame(Frame):# 继承Frame类
	def __init__(self,master):
		Frame.__init__(self,master)
		self.root = master #定义内部变量root
		self.createPage()

	def createPage(self):
		Label(self,text='开发者：梁凡\n班级：计算机161\n学号：2016123163', justify = LEFT, font = ('宋体', 15)).pack()
