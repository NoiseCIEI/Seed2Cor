#!/usr/bin/env python
"""
Wrapper for Seed2Cor.
"""
import os
from os import system
import sys
import logging.config
import textwrap

import yaml

from pymodule.timing import timing
from pymodule.sys_tool import setup_logging, mkdir


logger = logging.getLogger()
setup_logging('/projects/shzh3924/Src/PyModule/data/logging.yml')


class S2CWrapper(object):
    def __init__(self, fconfig="config.yml"):
        """
        :param fconfig: YAML file containing all parameters for Seed2Cor
        """
        self.fconfig = fconfig

        with open(self.fconfig, "r") as fcfg:
            self.cfg = yaml.safe_load(fcfg)

        cg = self.cfg
        self.lst_sta = os.path.join(cg['dir_seed'], 'station.lst')
        self.lst_seed = os.path.join(cg['dir_seed'], 'seed.lst')
        self.fpar = os.path.join(os.path.dirname(cg['dir_seed']), cg['fpar'])

    def _mklst(self, sh_sta='generate_stationlst.sh',
               sh_seed='generate_seedlst.sh'):
        """
        Make lists of stations & SEEDS.
        """
        dir_seed = self.cfg['dir_seed']

        if (os.path.exists(self.lst_sta)
            and os.path.exists(self.lst_seed)
            and self.cfg['fskip']['mklst']):
            logger.info('Skip making lists.')
            return

        for sh in [sh_sta, sh_seed]:
            system(f'{sh} {dir_seed}')

        return

    def _wfpar(self):
        """
        Write parameter file.
        """
        cg = self.cfg
        if os.path.exists(self.fpar) and cg['fskip']['wpar']:
            logger.info(f"Skip writing {cg['fpar']}")
            return

        ct = '\n'.join([
            f'{cg["rdseed"]}', f'{cg["evalresp"]}',
            self.lst_sta, self.lst_seed,
            f'{cg["cha"]}', f'{cg["sr"]}', f'{cg["gap_frac"]}', f'{cg["t1"]}',
            f'{cg["tlen"]}', f'{cg["Tmin"]}', f'{cg["Tmax"]}',
            f'{cg["fproc"]["tnorm"]}',
            f'{cg["eTmin"]}', f'{cg["eTmax"]}', f'{cg["ram_hlen"]}',
            f'{cg["sw_hlen"]}', f'{cg["fsmz"]}', f'{cg["fproc"]["tlen"]}',
            f'{cg["fproc"]["prc"]}', f'{cg["max_mem"]}',
            f'{cg["lag"]}', f'{cg["min_tlen"]}',
            f'{cg["fdel"]["sac"]}', f'{cg["fdel"]["amph"]}',
            f'{cg["fskip"]["ext_sac"]}', f'{cg["fskip"]["rm_resp"]}',
            f'{cg["fskip"]["norm"]}', f'{cg["fskip"]["cor"]}',
            f'{cg["fout"]}',
            ])

        with open(self.fpar, "w") as f:
            f.write(ct)

        return

    def _s2c(self, out='COR'):
        """
        Run Seed2Cor.
        """
        cg = self.cfg
        dir_out = os.path.join(os.path.dirname(cg['dir_seed']), out)
        mkdir(dir_out)
        os.chdir(dir_out)

        cmd = textwrap.dedent(
            f"""
            Seed2Cor {self.fpar} {cg['nthread']} <<- END
            Y
            END
            """
            )
        system(cmd)

        return

    def flow(self):
        """
        Workflow.
        """
        logger.debug('Make lists of stations & SEEDS.')
        self._mklst()
        logger.debug(f'Convert {self.fconfig} to {self.cfg["fpar"]}.')
        self._wfpar()
        logger.debug('Run Seed2Cor.')
        self._s2c()

        if self.cfg["fdel"]["par"]:
            logger.debug(f'Remove self.cfg["fpar"]')
            os.remove(self.fpar)

        return


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: python {sys.argv[0]} dirname param.yml")
    os.chdir(sys.argv[1])

    timing()
    wrapper = S2CWrapper(fconfig=sys.argv[2])
    wrapper.flow()
