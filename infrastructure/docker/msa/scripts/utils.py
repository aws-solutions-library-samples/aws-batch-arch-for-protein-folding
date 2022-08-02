# Copyright 2021 AlQuraishi Laboratory
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
from datetime import date


def add_data_args(parser: argparse.ArgumentParser):
    parser.add_argument(
        '--uniref90_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--mgnify_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--pdb70_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--uniclust30_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--bfd_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--jackhmmer_binary_path', type=str, default='/usr/bin/jackhmmer'
    )
    parser.add_argument(
        '--hhblits_binary_path', type=str, default='/usr/bin/hhblits'
    )
    parser.add_argument(
        '--hhsearch_binary_path', type=str, default='/usr/bin/hhsearch'
    )
    parser.add_argument(
        '--kalign_binary_path', type=str, default='/usr/bin/kalign'
    )
    parser.add_argument(
        '--max_template_date', type=str, 
        default=date.today().strftime("%Y-%m-%d"),
    )
    parser.add_argument(
        '--obsolete_pdbs_path', type=str, default=None
    )
    parser.add_argument(
        '--release_dates_path', type=str, default=None
    )
