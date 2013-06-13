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
\usepackage[table]{xcolor}

\usepackage{longtable}

\usepackage{pdflscape}

\begin{document}

\setlength\LTleft{0pt}
\setlength\LTright{0pt}
\setlength\tabcolsep{2pt}

{% for gene, variants, coverage in objects %}
    \begin{landscape}
    \section*{((( gene.name ))) gene}

    Report for the ((( gene.name ))) gene, accession number {\tt ((( gene.accession|escape_tex )))}, at position {\tt ((( gene.chrom ))):((( gene.start|int_add_commas )))--((( gene.end|int_add_commas )))}.

    {% for table in tables %}
        \subsection*{sss((( table.title|escape_tex )))}
        {\sffamily
            \begin{longtable}{@{\extracolsep{\fill}}((( table.column_alignments )))@{}}
            \toprule
                {% for column in table %}
                    ((( column.name|escape_tex ))){% if not loop.last %}&{% endif %}
                {% endfor %} \\
            \midrule
            \endhead
                {% for variant in variants() -%}
                    {% set call = variant.get_call(sample) %}
                    {% if call.is_variant and table.matcher(variant, call) -%}
                        {% for column in table -%}
                            {% if column.type == 'exon' -%}
                                E((( variant.exon )))
                            {% elif column.type == 'pos' %}
                                ((( variant.exon_offset )))
                            {% elif column.type == 'abspos' %}
                                \tt ((( variant.CHROM ))):((( variant.POS|int_add_commas )))
                            {% elif column.type == 'qual' %}
                            ((( call.data.GQ )))
                            {% elif column.type == 'type' %}
                            ((( variant.var_type )))
                            {% elif column.type == 'nc' %}
                                {% if variant.get_base_change(call)|join(", ")|length > 10 %}\footnotesize{% endif %}
                                {% if variant.get_base_change(call)|join(", ")|length > 15 %}\scriptsize{% endif %}
                                ((( variant.REF ))) $\rightarrow$ ((( variant.get_base_change(call)|join(", ") )))
                                    ({% if call.is_het %}het{% else %}hom{% endif %})
                            {% elif column.type == 'aac' %}
                                {% if variant.aa_change_texts %}
                                \tt ((( variant.aa_change_texts|join(", ")|escape_tex )))
                                {% else %}
                                    ---
                                {% endif %}
                            {% elif column.type == 'aaf' %}
                                ((( variant.aaf )))
                            {% elif column.type == 'webref' %}
                                {% if variant.ID -%}
                                    \href{http://ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=((( variant.ID )))}{((( variant.ID )))}
                                {% else %}
                                    ---
                                {%- endif %}
                            {% elif column.type == 'mutref' %}
                                {% if variant.get_mut_ref(call)|length > 20 %}\footnotesize{% endif %}
                                {% if variant.get_mut_ref(call)|length > 30 %}\scriptsize{% endif %}
                                ((( variant.get_mut_ref(call)|escape_tex )))
                            {% elif column.type == 'eff' %}
                                {% if not variant.effects -%}
                                    ---
                                {% else %}
                                    {% for effect in variant.effects -%}
                                        {% if effect.impact == 'HIGH' %}\color{red}{% endif %}
                                        ((( effect.effect_type_text ))){% if not loop.last %}, {% endif %}
                                    {%- endfor %}
                                {%- endif %}
                            {%- endif %}
                            {% if not loop.last %}&{% endif %}
                        {%- endfor %} \\
                    {%- endif %}
                {%- endfor %}
            \bottomrule
            \end{longtable}
        }
    {% endfor %}
    \end{landscape}
{% endfor %}
\end{document}
""")

    return template