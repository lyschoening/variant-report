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
\documentclass[a4paper,landscape]{report}
\usepackage{a4wide}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{microtype}

\usepackage{color}
\usepackage{tikz}
\usepackage{pgfplots}
\usetikzlibrary{datavisualization}
\usetikzlibrary{pgfplots.groupplots}

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

{% for gene, exons, max_pileup, coverage_table_rows in objects %}
    % \begin{landscape}
    \section*{((( gene.name ))) gene}

    Report for the ((( gene.name ))) gene, accession number {\tt ((( gene.accession|escape_tex )))}, at position {\tt ((( gene.chrom ))):((( gene.start|int_add_commas )))--((( gene.end|int_add_commas )))}.\\

    \begin{tikzpicture}

        \path (0cm,0cm) -- (\textwidth, 0cm);

        \begin{groupplot}[
            group style={
                columns=((( exons|length ))),
                vertical sep=0pt,
                ylabels at=edge left,
                xlabels at=edge bottom,
                yticklabels at=edge left,
                horizontal sep=0pt,
            },
            extra y ticks={20},
            extra y tick style={grid=major},
            width=1/(((exons|length )))*(\textwidth-2.0cm),
            height=3.0cm,
            ylabel=coverage,
            tickpos=left,
            xtick=\empty,
            ytick align=inside,
            xtick align=inside,
            scale only axis,
            ymin=0,
            ymax=((( max_pileup )))
        ]
        {% for exon in exons %}
        \nextgroupplot[xlabel=((( exon.name )))]
        \addplot [{% if exon.coding %}blue{% else %}black{% endif %}!80!black, fill={% if exon.coding %}blue{% else %}black{% endif %}, fill opacity=0.2] coordinates {{% for point in exon.points %}( ((( point.0 ))),((( point.1 ))) ){% endfor %}}
            |- (axis cs:0,0) -- cycle;
        {% endfor %}
        \end{groupplot}
    \end{tikzpicture}

    {\sffamily
        \begin{longtable}{@{\extracolsep{\fill}}llrrrrrrrrrrrrrrrrrrrr@{}}


        {% for indexes, exon_names in coverage_table_rows %}
                {% if loop.first %}
                    \toprule
		            & & \multicolumn{20}{c}{Exon} \\
                {% else %}
                    \midrule
                {% endif %}

                &
		        {% for name in exon_names %}
                    & ((( name )))
		        {% endfor %}\\

		        \midrule

                {% if loop.first %}
                    \endhead
                {% endif %}

                Coverage & avg.
                {% for index in indexes %}
                    & ((( '%d' | format(exons[index].mean) )))
                {% endfor %}\\
                 & min.
                {% for index in indexes %}
                    & ((( '%d' | format(exons[index].min) )))
                {% endfor %}\\
                 & max.
                {% for index in indexes %}
                    & ((( '%d' | format(exons[index].max) )))
                {% endfor %}\\
        {% endfor %}
        \bottomrule
		\end{longtable}
    }
    % \end{landscape}
{% endfor %}
\end{document}
""")

    return template