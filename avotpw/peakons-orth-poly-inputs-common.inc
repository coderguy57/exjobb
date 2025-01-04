% Common inputs for arXiv and Elsevier.

\usepackage{microtype}

\usepackage{amsmath}
\usepackage{amsthm}
\usepackage{amssymb}
\numberwithin{equation}{section}

\usepackage[dvipsnames]{xcolor}

\usepackage{tabularx}

\usepackage[showseconds=false]{datetime2}

\theoremstyle{plain}
\newtheorem{proposition}{Proposition}[section]
\newtheorem{corollary}[proposition]{Corollary}
\newtheorem{lemma}[proposition]{Lemma}
\newtheorem{theorem}[proposition]{Theorem}

\theoremstyle{definition}
\newtheorem{definition}[proposition]{Definition}
\newtheorem{example}[proposition]{Example}
\newtheorem{remark}[proposition]{Remark}
\newtheorem{assumption}[proposition]{Assumption}

% \usepackage{tikz}

% \usepackage[nottoc,notlot,notlof]{tocbibind} % To get "References" in the table of contents

\newcommand{\TODO}[1]{\textbf{\textup{\color{Magenta}[TODO:~{#1}]}}}
\newcommand{\HL}[1]{\textbf{\textup{\color{WildStrawberry}[HL:~{#1}]}}}
\newcommand{\JS}[1]{\textbf{\textup{\color{OrangeRed}[JS:~{#1}]}}}

\newcommand{\R}{\mathbf{R}}
\newcommand{\C}{\mathbf{C}}
\newcommand{\Z}{\mathbf{Z}}
\newcommand{\N}{\mathbf{N}}
\renewcommand{\epsilon}{\varepsilon}
\renewcommand{\phi}{\varphi}

\newcommand{\floor}[1]{\lfloor #1 \rfloor}
\newcommand{\abs}[1]{\left\lvert #1 \right\rvert}
\newcommand{\avg}[1]{\bigl\langle #1 \bigr\rangle}
\newcommand{\jump}[1]{\bigl[ #1 \bigr]}
\newcommand{\twin}{\widetilde}
\newcommand{\da}{d\alpha(x)}
\newcommand{\db}{d\beta(y)}
\newcommand{\dA}[1]{d\alpha^{#1}(x)}
\newcommand{\dB}[1]{d\beta^{#1}(y)}
\newcommand{\heineintegral}{\mathcal{J}}
\newcommand{\piecewiseclass}{PC^{\infty}}
\newcommand{\spaceDprime}{\mathcal{D}'(\R)}

% \newcommand{\circleS}{\mathbf{S}^1}
\newcommand{\circleS}{S^1}

\DeclareMathOperator{\sgn}{sgn}
\DeclareMathOperator*{\res}{res}

\DeclareMathOperator{\ad}{ad}
\DeclareMathOperator{\diff}{diff}
\DeclareMathOperator{\Diff}{Diff}

\hyphenation{self-adjoint anti-peakon anti-peakons}



%%% Local Variables:
%%% mode: latex
%%% TeX-master: "peakons-orth-poly-arXiv"
%%% End:
