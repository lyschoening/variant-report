from jinja2 import Environment
import re

__author__ = 'lyschoening'

LATEX_SUBS = (
    (re.compile(r'\\'), r'\\textbackslash'),
    (re.compile(r'([{}_#%&$])'), r'\\\1'),
    (re.compile(r'~'), r'\~{}'),
    (re.compile(r'\^'), r'\^{}'),
    (re.compile(r'"'), r"''"),
    (re.compile(r'\.\.\.+'), r'\\ldots'),
    )


def escape_tex(value):
    newval = value
    for pattern, replacement in LATEX_SUBS:
        newval = pattern.sub(replacement, newval)
    return newval

import locale
locale.setlocale(locale.LC_ALL, 'en_US')

def int_add_commas(number):
    return locale.format("%d", number, grouping=True)

def get_template():
    texenv = Environment()
    #    texenv.block_start_string = '((*'
    #    texenv.block_end_string = '*))'
    texenv.variable_start_string = '((('
    texenv.variable_end_string = ')))'
    #    texenv.comment_start_string = '((='
    #    texenv.comment_end_string = '=))'
    texenv.filters['escape_tex'] = escape_tex
    texenv.filters['int_add_commas'] = int_add_commas


    template = texenv.from_string(r"""
\documentclass[a4paper]{report}
\usepackage{a4wide}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{microtype}

\usepackage{color}

\usepackage[T1]{fontenc}
\usepackage{lmodern}

\usepackage{booktabs}
\heavyrulewidth=.13em
\lightrulewidth=.08em
\cmidrulewidth=.03em

\usepackage{hyperref}

\usepackage{longtable}

\usepackage{pdflscape}

\begin{document}

\setlength\LTleft{0pt}
\setlength\LTright{0pt}
\setlength\tabcolsep{2pt}

\title{((( sample_name|escape_tex )))}

{	\sf
    \begin{longtable}{@{\extracolsep{\fill}}llrrr@{\extracolsep{\fill}}}
    \toprule
    ID & Position & Min. & Avg. & Max. Coverage \\
    \endhead
    \toprule

{% for id contig, start, end, min, max, mean in objects %}
    \textbf{((( id|escape_tex )))} &
    \tt ((( gene.chrom ))):((( gene.start|int_add_commas )))--((( gene.end|int_add_commas ))) &
    ((( min ))) &
    ((( mean ))) &
    ((( max ))) \\
{% endfor %}

    \bottomrule
    \endfoot

    \end{longtable}
}

\end{document}
""")

    return template