#ifndef DATABASE_H_
#define DATABASE_H_

#include "lib.h"

struct Database {
  sqlite3 *db;
  Database() {
    db = NULL;
  }
  Database(string database_name) {
    string command = "sqlite3 " + database_name + " < create_tables.sql;";
    ASSERT(
        !system(command.c_str()) and sqlite3_open(database_name.c_str(), &db) == SQLITE_OK,
        "failed to open database"
        );
    // asynchronous to be more efficient
    sqlite3_exec(db, "PRAGMA synchronous = OFF", 0, 0, 0);
  }
  ~Database() {
    ASSERT(
        sqlite3_close(db) == SQLITE_OK,
        "failed to close database"
        );
  }
  // prepare for batch
  void open_transaction() {
    sqlite3_exec(db, "BEGIN TRANSACTION", 0, 0, 0);
  }
  void close_transaction() {
    sqlite3_exec(db, "END TRANSACTION", 0, 0, 0);
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
