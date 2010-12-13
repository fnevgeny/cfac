TOP  = ..

include $(TOP)/Make.lib

dbase_sql.o : schema.i

schema.i : schema.sql
	./c2cstr.sh < $? > $@
