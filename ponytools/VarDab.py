from minus80 import Freezable
from itertools import chain,repeat

import pandas as pd

class VarDab(Freezable):

    def __init__(self,name):
        super().__init__(name)
        self._initialize_tables()


    def genotypes(self,gformat='GT'):
        return pd.DataFrame(
            self._db.cursor().execute(
                'SELECT * FROM sample_genotypes WHERE format = ?',
                (gformat,)
            ).fetchall(),
            columns = ['chrom','pos','alt','ref','sample','format','value']
        )

    def add_VCF(self, filename):
        cur = self._db.cursor()
        cur.execute('BEGIN TRANSACTION')
        try:
            # Add the filename to the Database
            cur.execute('''
                INSERT OR IGNORE INTO files (filename) VALUES (?)
            ''',(filename,))
            file_id, = cur.execute('SELECT FILEID FROM files WHERE filename = ?;',(filename,)).fetchone()
            # Get a list of current variants and their ids
            cur_vars = { 
                (chrom,pos):(VID,id) \
                    for VID,chrom,pos,id \
                    in cur.execute(
                        'SELECT VARIANTID, chrom, pos, id FROM variants'
                    )
            }
            # Iterate over the file and build the pieces of the database
            with open(filename,'r') as IN:
                for line_id,line in enumerate(IN):
                    line = line.strip()
                    # Case 1: Header line
                    if line.startswith('##'):
                        items = self.parse_header(line)
                        for item in items:
                            cur.execute(
                                "INSERT INTO headers VALUES ({},{},?,?,?)".format(file_id,line_id),
                                item
                            )
                            if item[0] == 'FORMAT' and item[1] == 'ID': 
                                cur.execute('INSERT INTO format (format) VALUES (?)',(item[2],))
                            elif item[0] == 'INFO' and item[1] == 'ID':
                                cur.execute('INSERT INTO info (info) VALUES (?)',(item[2],))
                    # Case 2: Sample Line
                    elif line.startswith('#CHROM'):
                        samples = line.split('\t')[9:]
                        cur.executemany(
                            'INSERT INTO samples (FILEID,sample) VALUES ({},?)'.format(file_id),
                            [(x,) for x in samples ]
                        )
                        # Parse out the INFO and FMT fields
                        cur_info = {
                            key:ID for key,ID in cur.execute('SELECT info,INFOID FROM info')        
                        } 
                        cur_format = {
                            key:ID for key,ID in cur.execute('SELECT format,FMTID FROM format')        
                        } 
                        cur_samples = {
                            key:ID for key,ID in cur.execute(
                                '''SELECT sample,SAMPLEID from samples 
                                   WHERE FILEID = "{}" 
                                   ORDER BY SAMPLEID'''.format(file_id)
                            )        
                        } 
                    # Case 3: Genotypes
                    else:
                        chrom,pos,id,ref,alt,qual,fltr,info,fmt,*genos = line.split('\t')
                        # Insert the variant
                        if (chrom,pos) not in cur_vars:
                            cur.execute('''
                                INSERT INTO variants (chrom,pos,id,ref,alt) VALUES (?,?,?,?,?)
                            ''',(chrom,pos,id,ref,alt))
                            var_id = self._db.last_insert_rowid()
                        else:
                            var_id = cur_vars[(chrom,pos)]
                        # Insert the observed QUAL score
                        cur.execute('''
                            INSERT INTO variant_qual (FILEID,VARIANTID,qual,filter) VALUES (?,?,?,?)
                        ''',(file_id,var_id,qual,fltr))
                        # Insert the info fields
                        info = [field.split('=') for field in info.split(';')]
                        fmt = [cur_format[x] for x in fmt.split(':')] 
                        # Insert the genotypes 
                        genos = [field.split(':') for field in genos]
                        records = []
                        for geno,sample in zip(genos,cur_samples.values()):
                            for i,g in enumerate(geno):
                                records.append((
                                    file_id,
                                    var_id,
                                    sample,
                                    fmt[i],
                                    g
                                ))
                        cur.executemany('''
                            INSERT INTO genotypes VALUES (?,?,?,?,?) 
                        ''',records)
            cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK;')
            raise e

    def _drop_tables(self):
        tables = [x[0] for x in self._db.cursor().execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall() ]
        for table in tables:
            if table != 'sqlite_sequence':
                self._db.cursor().execute('''
                    DROP TABLE IF EXISTS {};
                '''.format(table))
        views = [x[0] for x in self._db.cursor().execute(
            "SELECT name FROM sqlite_master WHERE type='view'"
        ).fetchall() ]
        for table in tables:
            if table != 'sqlite_sequence':
                self._db.cursor().execute('''
                    DROP VIEW IF EXISTS {};
                '''.format(table))



    def _initialize_tables(self):
        cur = self._db.cursor()
        # Files -- Meta Data
        cur.execute('''
        CREATE TABLE IF NOT EXISTS files (
            FILEID INTEGER PRIMARY KEY AUTOINCREMENT,
            filename TEXT NOT NULL UNIQUE,
            added datetime DEFAULT CURRENT_TIMESTAMP
        );
        ''')
        # Headers
        cur.execute('''
        CREATE TABLE IF NOT EXISTS headers (
            FILEID INT,  -- Maps to files
            LINEID INT,  -- Records the lin
            tag TEXT,
            key TEXT DEFAULT NULL,
            val TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID)
        );
        ''')
        # Samples
        cur.execute('''
        CREATE TABLE IF NOT EXISTS samples (
            FILEID INTEGER,  -- Maps to files
            SAMPLEID INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID)
        );
        ''')
        # Variants
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variants (
            VARIANTID INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT,
            pos INT,
            id TEXT,
            ref TEXT,
            alt TEXT
        );
        ''')
        # Variant Quality
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variant_qual (
            FILEID INTEGER,
            VARIANTID INTEGER,
            qual REAL,
            filter TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID)
        );
        ''')
        # Variant Info
        cur.execute('''
        CREATE TABLE IF NOT EXISTS info (
            INFOID INTEGER PRIMARY KEY AUTOINCREMENT,
            info TEXT
        );
        ''')
        # Variant Info
        cur.execute('''
        CREATE TABLE IF NOT EXISTS variant_info(
            FILEID INTEGER,
            INFOID INTEGER,
            VARIANTID INTEGER,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(INFOID) REFERENCES info(INFOID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID)
        );
        ''')
        # FORMAT
        cur.execute('''
        CREATE TABLE IF NOT EXISTS format (
            FMTID INTEGER PRIMARY KEY AUTOINCREMENT,
            format TEXT UNIQUE
        );
        ''')
        # Genotypes
        cur.execute('''
        CREATE TABLE IF NOT EXISTS genotypes (
            FILEID INT,
            VARIANTID INT,
            SAMPLEID INT,
            FMTID INT,
            value TEXT,
            FOREIGN KEY(FILEID) REFERENCES files(FILEID),
            FOREIGN KEY(VARIANTID) REFERENCES variants(VARIANTID),
            FOREIGN KEY(SAMPLEID) REFERENCES samples(SAMPLEID),
            FOREIGN KEY(FMTID) REFERENCES format(FMTID)
        ); 
        ''')

        cur.execute(''' 
            CREATE VIEW IF NOT EXISTS sample_genotypes AS
            SELECT 
             chrom,pos,alt,ref,
             sample, format, value 
            FROM genotypes 
            INNER JOIN variants on genotypes.VARIANTID = variants.VARIANTID
            INNER JOIN samples on genotypes.SAMPLEID = samples.SAMPLEID
            INNER JOIN format on genotypes.FMTID = format.FMTID
        ''')

    # -----------------------------------------------------------------
    # Static Methods
    # -----------------------------------------------------------------

    def parse_header(self,header_string):
        '''
        Parses a header line into a list of tuples: (tag, key, val)
        Example
        -------
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        gets parsed into:
            [('FORMAT', 'ID', 'GQ'),
             ('FORMAT', 'Number', '1'),
             ('FORMAT', 'Type', 'Integer'),
             ('FORMAT', 'Description', '"Genotype Quality"')]

        '''
        header_string = header_string.strip('#')
        tag = None
        keys = []
        vals = []
        current = ''
        in_str = False
        str_delim = None
        for i,char in enumerate(header_string):
            if in_str:
                if (char == '"' or char == '"') and header_string[i-1] != '\\': 
                    current += char
                    vals.append(current)
                    current = ''
                    in_str = False
                else:
                    current += char
            elif char == '"' or char == "'":
                current += char
                in_str = True
            elif char == '=':
                if tag == None:
                    tag = current
                else:
                    if len(keys) == len(vals):
                        keys.append(current)
                    else:
                        raise ValueError('fuck')
                current = ''
            elif char == ',' or char == '>' or char == '<':
                if in_str:
                    current += char
                elif len(vals) == len(keys)-1:
                    vals.append(current)
                    current = ''
            else:
                current += char 
        # Clean up the loose ends
        if current != '':
            vals.append(current)
        if len(vals) == 1:
            keys.append(None)
        return [ (tag,x,y) for x,y in zip(keys,vals)]

