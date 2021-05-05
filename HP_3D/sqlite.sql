-- 创建表result, 用于存放程序结果
-- id, 整型, 主键, 自动递增
-- fitness, 整型, 适应度
-- seq, TEXT, 摆放序列
create table result(
	id integer primary key autoincrement,
	fitness int not null,
	seq text not null
);

-- 创建表result_log, 用于记录日志信息
-- log_id, 来自表result中的id, 要与result.id数据类型保持一致
-- operate, 记录insert, update, delete操作
-- entry_datetime, 记录被创建时的时间
create table result_log(
	log_id integer not null,
	operate varchar(10) not null,
	entry_datetime text not null
);

-- 在result表上创建一个触发器result_insert_trigger
-- 当result表中有数据插入时，在result_log表中添加相应的日志
create trigger result_insert_trigger after insert
on result
begin
	insert into result_log (log_id, operate, entry_datetime) values (new.id, 'insert', datetime('now'));
end;

-- 在result表上创建一个触发器result_update_trigger
-- 当result表中有数据更新时，在result_log表中更新相应的日志
create trigger result_update_trigger after update
on result
begin
	update result_log set operate = 'update', entry_datetime = datetime('now') where log_id = new.id;
end;

-- 在result表上创建一个触发器result_delete_trigger
-- 当result表中有数据删除时，在result_log表中更新相应的日志
create trigger result_delete_trigger before delete
on result
begin
	update result_log set operate = 'delete', entry_datetime = datetime('now') where log_id = old.id;
end;

-- 从result表创建视图: 获取最优值及其序列
create view result_best_fitness_view as
select min(fitness) as best_fitness, seq 
from result;

-- 从result表创建视图: 获取最差值及其序列
create view result_worst_fitness_view as
select max(fitness) as worst_fitness, seq
from result;

-- 创建一个基于表result的列fitness上的索引
create index result_fitness_index 
on result(fitness);

-- 创建表user, 用于存放用户登录信息
-- uid, 整型, 主键, 自动递增
-- uname, TEXT, 账号
-- upassword, TEXT, 密码
create table user(
	uid integer primary key autoincrement,
	uname text not null,
	upassword text not null
);
