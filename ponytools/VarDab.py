from minus80 import Freezable
import re

class VarDab(Freezable):

    def __init__(self,name):
        super().__init__(name)
        self._initialize_tables()


    def add_VCF(self, filename):
        cur = self._db.cursor()
        try:
            with open(filename,'r') as IN:
                cur.execute('BEGIN TRANSACTION')
                cur.execute('''
                    INSERT INTO files (filename) VALUES (?)
                ''',(filename,))
                file_id, = cur.execute('SELECT FILEID FROM files WHERE filename = ?;',(filename,)).fetchone()
                for line in IN:
                    if line.startswith('##'):
                        import ipdb; ipdb.set_trace()
                        continue
                        tag,val = line.strip().replace('##','').split('=',maxsplit=1)
                        if val.startswith('<'):
                            for item in val.replace('>','').split(','):
                                key,val = item.split('=')
                                cur.execute('''
                                    INSERT OR REPLACE INTO headers (FILEID, tag, key, val) 
                                    VALUES (?,?,?,?)
                                ''',(file_id,tag,key,val))
                        else:
                            cur.execute('''
                                INSERT OR REPLACE INTO headers (FILEID,tag,val)
                                VALUES (?,?,?)
                            ''',(file_id,tag,val))
                    else:
                        continue
                cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK;')
            raise e

    def _drop_tables(self):
        self._db.cursor().execute('''
            DROP TABLE files;
            DROP TABLE headers;
            DROP TABLE variants;
        ''')


    def _initialize_tables(self):
        cur = self._db.cursor()
        cur.execute('''
        CREATE TABLE IF NOT EXISTS files (
            FILEID INTEGER PRIMARY KEY AUTOINCREMENT,
            filename TEXT NOT NULL UNIQUE,
            added datetime DEFAULT CURRENT_TIMESTAMP
        );
        CREATE TABLE IF NOT EXISTS headers (
            FILEID INTEGER,  -- Maps to files
            tag TEXT,
            key TEXT DEFAULT NULL,
            val TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID)
        );
        CREATE TABLE IF NOT EXISTS variants (
            VARIANTID INTEGER PRIMARY KEY,
            chrom TEXT,
            pos INT,
            id TEXT,
            ref TEXT,
            alt TEXT,
            qual REAL, 
            filter TEXT
            )
        ''')
