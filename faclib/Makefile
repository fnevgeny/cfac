TOP  = ..

include $(TOP)/Make.lib

dbase_sql.o : schema.i

schema.i : schema.sql
	./sql2cstr.sh < $? > $@
