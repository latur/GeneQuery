# GeneQuery

## Command line version [gq.py]

> Requirements: python3, genequery

### Install:

```bash
wget 'https://...'
unzip default.gqdb

wget 'https://...'
```

### Usage:

```
Usage:
  python3 gq.py -g <genes> -s <sp>:<sp> [-db <path>]

Arguments:
  -g  <genes>    Gene names, separated by space or comma
  -s  <sp>:<sp>  Query species shortcode : Database species shortcode
  -db <path>     Qustom database path.  Default: default.gqdb

Examples:
  python3 gq.py -g "15368 93692 ... 104263 99929"
  python3 gq.py -g 15368,93692,...,104263,99929
  python3 gq.py -s mm:rn -g 15368,93692,...,104263,99929
  python3 gq.py -s gp:gp -db custom.gqdb -g 93692,15368,...,104263
```

## Web Server version [gqserver.py]

> Requirements: python3, genequery, Flask

### Install:

```bash
wget 'https://...'
unzip default.gqdb

git clone https://github.com/latur/GeneQuery
cd GeneQuery
```

### Run web-server:

```
Usage:
  python3 gqserver.py [-p <port>] [-db <path>]

Arguments:
  -p  <port>  Server port.           Default: 9225
  -db <path>  Qustom database path.  Default: default.gqdb

Examples:
  python3 gqserver.py
  python3 gqserver.py -db custom.gqdb -p 80
  python3 gqserver.py -h
```

Open in your browser http://0.0.0.0:9225/


## How to create your own GeneQuery Database
