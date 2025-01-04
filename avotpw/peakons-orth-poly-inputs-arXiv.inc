% Inputs for arXiv version

\input{peakons-orth-poly-inputs-common.inc}

\usepackage{newtxtext,newtxmath}

\usepackage{graphicx}

\usepackage[margin=10pt,font=small,labelfont=bf,labelsep=period]{caption}

\usepackage[nottoc]{tocbibind}

\usepackage{titlesec}
\titleformat*{\section}{\large\bfseries}
\titleformat*{\subsection}{\normalsize\bfseries}

\usepackage[twocolumn, columnsep=7mm, left=19mm, right=14mm, top=25mm, bottom=30mm]{geometry}

\usepackage[colorlinks,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}


\newcommand{\makefrontmatter}{%
  \hypersetup{%
    pdfauthor={Hans Lundmark, Jacek Szmigielski},
    pdftitle={\ourtitle},
  }

  \begin{center}
    {\LARGE \ourtitle}

    \vspace{8mm}
  
    \begin{minipage}{13.5cm}
      \begin{minipage}[t]{4.5cm}
        \large
        \kern0pt
        Hans Lundmark
      \end{minipage}
      \begin{minipage}[t]{9cm}
        \kern0pt
        Department of Mathematics \\
        Linköping University \\
        SE-581\,83 Linköping, Sweden \\
        hans.lundmark@liu.se 
      \end{minipage}
      
      \vspace{5mm}
      
      \begin{minipage}[t]{4.5cm}
        \large
        \kern0pt
        Jacek Szmigielski
      \end{minipage}
      \begin{minipage}[t]{9cm}
        \kern0pt
        Department of Mathematics and Statistics \& \\
        \hspace*{1em}Centre for Quantum Topology and Its Applications (quanTA) \\
        University of Saskatchewan \\
        106 Wiggins Road, Saskatoon, Saskatchewan, S7N\,5E6, Canada \\
        szmigiel@math.usask.ca
      \end{minipage}
    \end{minipage}

    \vspace{6mm}

    {\normalsize \ourdate}
  \end{center}

  % \vspace{1mm}

  \begin{quote}
    \small \textbf{Abstract.} \ourabstract{}
  \end{quote}
  
  \vspace{4mm}
  
  % Don't bother with keywords.
}

%%% Local Variables:
%%% mode: latex
%%% TeX-master: "peakons-orth-poly-arXiv"
%%% End:
