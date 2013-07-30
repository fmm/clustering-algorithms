#ifndef DATABASE_H_
#define DATABASE_H_

#include "lib.h"
#include <sqlite3.h>

struct Database {
	sqlite3 *db;
	Database() {
	  db = NULL;
	}
	Database(string database_name) {
		ASSERT(
			sqlite3_open(database_name.c_str(), &db) == SQLITE_OK,
			"failed to open database"
		);
	}
	~Database() {
		ASSERT(
		    sqlite3_close(db) == SQLITE_OK,
				"failed to close database"
		);
	}
	// the reference of the database is lost when you use the assignment operator
	Database& operator=(Database rhs) {
	  this->db = rhs.db, rhs.db = NULL;
	  return *this;
	}
	vector< vector<string> > execute(string sql) {
		ASSERT(db, "database is not opened");
		vector< vector<string> > result;
		sqlite3_stmt *statement;
		ASSERT(
				sqlite3_prepare_v2(db, sql.c_str(), -1, &statement, 0) == SQLITE_OK,
				"failed to execute sql=\"" + sql + "\""
		);
		for(int columns = sqlite3_column_count(statement); sqlite3_step(statement) == SQLITE_ROW;) {
			vector<string> values(columns, "?");
			for(int column = 0; column < columns; ++ column) {
				if(char *ptr = (char*)sqlite3_column_text(statement, column)) {
					values[column] = string(ptr);
				}
			}
			result.push_back(values);
		}
		ASSERT(
			sqlite3_finalize(statement) == SQLITE_OK,
			"failed to finalize statement"
		);
		return result;
	}
};
#endif
