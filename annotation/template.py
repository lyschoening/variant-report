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

    \subsection*{Unique, Moderate- \& High-impact Variants}
    {\sffamily
        \begin{longtable}{@{\extracolsep{\fill}}llrclrrlr@{}}
        \toprule
            Location & Position & Quality & Type & Nucleotide change & AAF & Web Ref (refSNP) & Mut Ref & Mutation Effect  \\
        \midrule
        \endhead
            {% for variant in variants() -%}
                {% set call = variant.get_call(sample) %}
                {% if call.is_variant and ((call.data.GQ > 9 and ('MODERATE' in call.impacts or variant.num_het + variant.num_hom_alt == 1)) or 'HIGH' in call.impacts) -%}
                    E((( variant.exon ))) &
                    ((( variant.exon_offset ))) &
                    ((( call.data.GQ ))) &
                    ((( variant.var_type ))) &
                        {% if variant.get_aa_change(call)|join(", ")|length > 10 %}\footnotesize{% endif %}
                        {% if variant.get_aa_change(call)|join(", ")|length > 15 %}\scriptsize{% endif %}
                        ((( variant.REF ))) $\rightarrow$ ((( variant.get_aa_change(call)|join(", ") )))
                            ({% if call.is_het %}het{% else %}hom{% endif %})
                     &
                    ((( variant.aaf ))) &
                    {% if variant.ID -%}
                        \href{http://ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=((( variant.ID )))}{((( variant.ID )))}
                    {% else %}
                        ---
                    {%- endif %} &
                    {% if variant.get_mut_ref(call)|length > 20 %}\footnotesize{% endif %}
                    {% if variant.get_mut_ref(call)|length > 30 %}\scriptsize{% endif %}
                    ((( variant.get_mut_ref(call)|escape_tex )))  &
                    {% if not variant.effects %}---{% endif %}
                    {% for effect in variant.effects -%}
                        {% if effect.impact == 'HIGH' %}\color{red}{% endif %}
                        ((( effect.effect_type_text ))){% if not loop.last %}, {% endif %}
                    {%- endfor %}\\
                {%- endif %}
            {%- endfor %}
        \bottomrule
        \end{longtable}
    }

    \subsection*{All Variants}
    {\sffamily
        \begin{longtable}{@{\extracolsep{\fill}}llrclrlr@{}}
        \toprule
            Location & Position & Quality & Type & Nucleotide change & Web Ref (refSNP) & Mut Ref & Mutation Effect  \\
        \midrule
        \endhead
            {% for variant in variants() -%}
                {% set call = variant.get_call(sample) %}
                {% if call.is_variant and (call.data.GQ > 9 or 'HIGH' in call.impacts) -%}
                    E((( variant.exon ))) &
                    ((( variant.exon_offset ))) &
                    ((( call.data.GQ ))) &
                    ((( variant.var_type ))) &
                        {% if variant.get_aa_change(call)|join(", ")|length > 10 %}\footnotesize{% endif %}
                        {% if variant.get_aa_change(call)|join(", ")|length > 15 %}\scriptsize{% endif %}
                        ((( variant.REF ))) $\rightarrow$ ((( variant.get_aa_change(call)|join(", ") )))
                            ({% if call.is_het %}het{% else %}hom{% endif %})
                     &
                    {% if variant.ID -%}
                        \href{http://ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=((( variant.ID )))}{((( variant.ID )))}
                    {% else %}
                        ---
                    {%- endif %} &
                    {% if variant.get_mut_ref(call)|length > 20 %}\footnotesize{% endif %}
                    {% if variant.get_mut_ref(call)|length > 30 %}\scriptsize{% endif %}
                    ((( variant.get_mut_ref(call)|escape_tex )))  &
                    {% if not variant.effects %}---{% endif %}
                    {% for effect in variant.effects -%}
                        {% if effect.impact == 'HIGH' %}\color{red}{% endif %}
                        ((( effect.effect_type_text ))){% if not loop.last %}, {% endif %}
                    {%- endfor %}\\
                {%- endif %}
            {%- endfor %}
        \bottomrule
        \end{longtable}
    }
    \end{landscape}
{% endfor %}
\end{document}
""")

    return template