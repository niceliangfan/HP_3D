from tkinter import *
from time import *
from tkinter.messagebox import *
from mainpage import *
from sqlite3 import *

class LoginPage(object):
	def __init__(self,master=None):
		self.root = master #定义内部变量root
		self.root.geometry('%dx%d' % (890,525)) #设置窗口大小
		self.username = StringVar()
		self.password = StringVar()
		self.nowtime = StringVar()
		self.createPage()

	def createPage(self):
		self.page = Frame(self.root) #创建Frame
		self.page.pack()
		Label(self.page, text = '欢迎进入系统', font = ('宋体',20)).grid(row=0, column = 1, stick=W)
		self.timelabel = Label(self.page, textvariable = self.nowtime, fg = 'blue', font = ('黑体', 12)).grid(row = 1, column = 2, stick = E)
		Label(self.page,text = '账号: ').grid(row=3,stick=W)
		Entry(self.page,textvariable=self.username).grid(row=3,column=1,stick=E)
		Label(self.page,text = '密码: ').grid(row=4,stick=W)
		Entry(self.page,textvariable=self.password,show='*').grid(row=4,column=1,stick=E)
		Label(self.page).grid(row=5,stick=W)
		Button(self.page,text='登录',command=self.loginCheck).grid(row=6,column=0,stick=W)
		Button(self.page,text='重置',command=self.reset).grid(row=6, column=1,stick=E)

		self.getTime()


	def loginCheck(self):
		flag = False
		name = self.username.get()
		secret = self.password.get()
		conn = connect('data.db')
		cursor = conn.cursor()
		cursor.execute('select uname, upassword from user')
		values = cursor.fetchall()
		cursor.close()
		conn.close()
		for value in values:
			if name == str(value[0]) and secret == str(value[1]):
				flag = True
				break
		if flag:
			if name == 'liangfan':
				self.page.destroy()
				AdminPage(self.root)
			else:
				self.page.destroy()
				UserPage(self.root)
		else:
			showinfo(title='登录失败',message='账号或密码错误！')

	def reset(self):
		self.username.set('')
		self.password.set('')

	def getTime(self):
		self.nowtime.set(strftime('%H:%M:%S'))
		self.page.after(1000, self.getTime)
