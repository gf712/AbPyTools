import configparser
from abpytools.home import Home

config = configparser.ConfigParser()
config.read(f"{Home().homedir}/config.ini")

HAS_PROTO = True if config["PROTOBUF"]['protobuf'] == '1' else False
