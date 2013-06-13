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

\usepackage{xcolor}
\usepackage{colortbl}

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
\maketitle

\begin{landscape}

{	\sf
    \begin{longtable}{@{\extracolsep{\fill}}lrrrr@{\extracolsep{\fill}}}
    \toprule
    ID & Position & Min. & Avg. & Max. Coverage \\
    \toprule
    \endhead

{% for name, chrom, start, end, min_coverage, max_coverage, mean_coverage in objects %}
    {% if min_coverage < minimum_coverage %}\color{red}{% endif %} \textbf{((( name|escape_tex )))} &
    \tt ((( chrom ))):((( start|int_add_commas )))--((( end|int_add_commas ))) &
    {% if min_coverage < minimum_coverage %}\cellcolor{red!25}{% endif %} ((( '%d' | format(min_coverage) ))) &
    {% if mean_coverage < minimum_coverage %}\cellcolor{red!25}{% endif %} ((( '%3.2f' | format(mean_coverage) ))) &
    {% if max_coverage < minimum_coverage %}\cellcolor{red!25}{% endif %} ((( '%d'| format(max_coverage) ))) \\
{% endfor %}

    \bottomrule
    \end{longtable}
}
\end{landscape}

\end{document}
""")

    return template