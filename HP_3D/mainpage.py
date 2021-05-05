from tkinter import *
from view import * #菜单栏对应的各个子页面

class AdminPage(object):
	def __init__(self,master=None):
		self.root = master #定义内部变量root
		self.root.geometry('%dx%d' % (890,525)) #设置窗口大小
		self.createPage()
     
	def createPage(self):
		self.operatePage = OperateFrame(self.root) # 创建不同Frame
		self.userPage = UserFrame(self.root)
		self.swInfoPage = SoftWareInfoFrame(self.root)
		self.dlInfoPage = DeveloperInfoFrame(self.root)
		self.operatePage.pack() #默认显示主界面

		mainmenu = Menu(self.root)

		interfacemenu = Menu(mainmenu)
		mainmenu.add_cascade(label = '界面', menu = interfacemenu)
		interfacemenu.add_command(label='操作',command = self.operate)
		interfacemenu.add_separator()
		interfacemenu.add_command(label = '退出', command = self.root.destroy)

		usermenu = Menu(mainmenu)
		mainmenu.add_cascade(label = '用户', menu = usermenu)
		usermenu.add_command(label = '用户管理', command = self.manageUser)

		aboutmenu = Menu(mainmenu)
		mainmenu.add_cascade(label = '关于', menu = aboutmenu)
		aboutmenu.add_command(label = '软件信息', command = self.showInfo)
		aboutmenu.add_command(label = '开发者信息', command = self.showDeveloper)
		self.root['menu'] = mainmenu # 设置菜单栏

	def operate(self):
		self.operatePage.pack()
		self.userPage.pack_forget()
		self.swInfoPage.pack_forget()
		self.dlInfoPage.pack_forget()

	def manageUser(self):
		self.operatePage.pack_forget()
		self.userPage.pack()
		self.swInfoPage.pack_forget()
		self.dlInfoPage.pack_forget()

	def showInfo(self):
		self.operatePage.pack_forget()
		self.userPage.pack_forget()
		self.swInfoPage.pack()
		self.dlInfoPage.pack_forget()

	def showDeveloper(self):
		self.operatePage.pack_forget()
		self.userPage.pack_forget()
		self.swInfoPage.pack_forget()
		self.dlInfoPage.pack()


class UserPage(object):
	def __init__(self,master=None):
		self.root = master #定义内部变量root
		self.root.geometry('%dx%d' % (890,525)) #设置窗口大小
		self.createPage()
     
	def createPage(self):
		self.operatePage = OperateFrame(self.root) # 创建不同Frame
		self.swInfoPage = SoftWareInfoFrame(self.root)
		self.dlInfoPage = DeveloperInfoFrame(self.root)
		self.operatePage.pack() #默认显示主界面

		mainmenu = Menu(self.root)

		interfacemenu = Menu(mainmenu)
		mainmenu.add_cascade(label = '界面', menu = interfacemenu)
		interfacemenu.add_command(label='操作',command = self.operate)
		interfacemenu.add_separator()
		interfacemenu.add_command(label = '退出', command = self.root.destroy)

		aboutmenu = Menu(mainmenu)
		mainmenu.add_cascade(label = '关于', menu = aboutmenu)
		aboutmenu.add_command(label = '软件信息', command = self.showInfo)
		aboutmenu.add_command(label = '开发者信息', command = self.showDeveloper)
		self.root['menu'] = mainmenu # 设置菜单栏

	def operate(self):
		self.operatePage.pack()
		self.swInfoPage.pack_forget()
		self.dlInfoPage.pack_forget()

	def showInfo(self):
		self.operatePage.pack_forget()
		self.swInfoPage.pack()
		self.dlInfoPage.pack_forget()

	def showDeveloper(self):
		self.operatePage.pack_forget()
		self.swInfoPage.pack_forget()
		self.dlInfoPage.pack()
