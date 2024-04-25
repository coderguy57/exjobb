%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Template for bachelor's and master's theses at the
%     Department of Mathematics at Linköping University.
%         version 0.9 (2016-05-11)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This template was made by Jonas Granholm
%     (jonas.b.granholm@gmail.com)
%
% It is a nearly complete rewrite of the previous
% template which was compiled by Per Erik Strandberg



% If you don't know what you're doing, don't edit this file.
% Create a file with the .tex extension containing the
% outline below, and start from there.
%
%     \documentclass[<language>,<level>]{liumaiex}
%     
%     <preamble>
%     
%     \begin{document}
%     
%     \frontmatter
%     
%     <information about thesis>
%     
%     \maketitle
%     
%     <introduction>
%     
%     \mainmatter
%     
%     <main text>
%     
%     \end{document}
%
% <language> should be replaced by either english or swedish.
% <level> should be replaced by either bachelor or master,
%   or not used at all.
% <preamble> is the place to load any packages you want to use
%   and perhaps define useful macros.
% <information about thesis> should look something like this:
%
%     \title{<the title of your thesis>}
%     \author{<your full name(s)>}
%     \shortauthor{<your last name(s)>}
%     \supervisor{<your supervisors name>}{<his/her department>}
%     \examiner{<your examiners name>}{<his/her department>}
%     \regnumber{<the ISRN given by administrators later>}
%     \publishmonth{<the month the thesis is published>}
%     \publishyear{<the year the thesis is published>}
%     \enkeywords{<the english keywords of your thesis, comma separated>}
%     \svkeywords{<the swedish keywords of your thesis, comma separated>}
%     \publishurl{<url to the published thesis (you'll get it later)>}
%
%   The command \supervisor can be used multiple times.
%   There are a few additional commands that may be used:
%
%     \department{<department, if not MAI>}
%     \city{<the city where the thesis is published, if not Linköping>}
%     \level{<usually G2 or A, if level is not specified earlier>}
%     \credits{<usually 16 hp or 30 hp, if level is not specified earlier>}
%     \pdftitle{<the title of the thesis>}
%     \pdfauthor{<the author(s) of the thesis>}
%     \pdfkeywords{<the keywords of the thesis>}
%     \pdfsubject{<the subject of the thesis>}
%
%   The information in the last four commands will end up in the
%   pdf metadata. If these commands are not used, then the title,
%   author and keywords info from the previous commands will be
%   used instead, but the \pdf... variants can be used if there
%   are problems with messy formatting etc.
% <introduction> should contain things like your abstract and
%   acknowledgement, and it should end with \tableofcontents
%   and perhaps \listoffigures and \listoftables.
%   The keywords and url are usually added after the abstract,
%   using the commands \placeenkeywords, \placesvkeywords,
%   \placeenurl and \placesvurl.
% <main text> is the rest of the thesis.
%
% To compile your thesis you need this file and either liulogo_en.eps
% or liulogo_sv.eps in the same directory as your .tex file.



% Copyright (c) 2016, Jonas Granholm
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation
% files (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
% THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


% This is the so called MIT license, and it means that you can do
% pretty much anything you want with the code, except for claiming
% that you wrote it or blaiming me if anything goes wrong.
%
% However, if you distribute a modified version of this template
% it is a good idea to clearly note that it differs from the original,
% preferably by changing the name of the class. That way it's less likely
% that people run into problems if they want to use the original template
% in the future.
%
% Note that you do not need to include this copyright text in documents
% compiled with the template. That is only necessary if you make a new
% LaTeX class based on this one.







%---------------------------------------
%	Identification
%---------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{liumaiex}[2016/05/11 Exjobbsmall för Matematiska institutionen vid Linköpings universitet v0.9]
% Matches perfectly with line length in tex log file, so change with caution


%---------------------------------------
%	Options
%---------------------------------------

% --- Language ---
\newif\ifliumaiex@langswe
% English is default since \newif sets new bools to false
% If both options are selected, a warning is given and English is used
\DeclareOption{swedish}{\liumaiex@langswetrue}
\DeclareOption{english}
	{\ifliumaiex@langswe
		\ClassWarningNoLine{liumaiex}
			{Both 'swedish' and 'english' options used, defaulting to english}
	\fi
	\liumaiex@langswefalse}

% --- Level ---
% If both options are selected, a warning is given and bachelor is used
\DeclareOption{bachelor}
	{\newcommand{\liumaiex@level}{G2}
	\newcommand{\liumaiex@credits}{16 hp}}
\DeclareOption{master}
	{\ifx\liumaiex@level\undefined
		\newcommand{\liumaiex@level}{A}
		\newcommand{\liumaiex@credits}{30 hp}
	\else
		\ClassWarningNoLine{liumaiex}{Both 'bachelor' and 'master' options used, defaulting to bachelor}
	\fi}

% --- Hyperref ---
\newif\ifliumaiex@hyperref
\liumaiex@hyperreftrue
\DeclareOption{nohyperref}{\liumaiex@hyperreffalse}

% --- Paper size ---
\newif\ifliumaiex@four
\DeclareOption{a4paper}{\liumaiex@fourtrue}


% --- Executing options ---
\ProcessOptions\relax


%---------------------------------------
%	Packages
%---------------------------------------

\LoadClass{book}

\RequirePackage{iftex}
\newif\ifliumaiex@xeorlua
\ifXeTeX\liumaiex@xeorluatrue\fi
\ifLuaTeX\liumaiex@xeorluatrue\fi
\ifliumaiex@xeorlua
	\RequirePackage{polyglossia}
	\RequirePackage{fontspec}
\else
	\RequirePackage[T1]{fontenc}
	\RequirePackage[swedish,english]{babel}
	\RequirePackage[utf8]{inputenc}
\fi
\RequirePackage{graphics}
\RequirePackage{epstopdf}
\RequirePackage{tabularx}
\RequirePackage{fancyhdr}
\ifliumaiex@hyperref
	\RequirePackage{hyperref}
\else
	% If hyperref is not loaded:
	\AtBeginDocument{\providecommand{\url}{\texttt}}
	% Not the best replacement, but it should work
	% for the link in \publishurl as well as for
	% the links on the copyright page.
\fi
\RequirePackage{geometry}
\ifliumaiex@four
	\geometry{textwidth=345pt,textheight=490pt,hcentering}
\else
	\geometry{paperwidth=165mm,paperheight=240mm,
		textwidth=345pt,textheight=490pt,top=3cm,inner=2cm}
\fi

%---------------------------------------
%	Info and labels
%---------------------------------------

% Info, will be set by author
\newcommand{\liumaiex@title}{}
\newcommand{\liumaiex@author}{}
\newcommand{\liumaiex@shortauthor}{}
\newcommand{\liumaiex@publishmonth}{}
\newcommand{\liumaiex@publishyear}{}
\newcommand{\liumaiex@city}{Linköping}
\newcommand{\liumaiex@department}{\mai}
\newcommand{\liumaiex@examiner}{{}{}}
\providecommand{\liumaiex@level}{} % If not set in options
\providecommand{\liumaiex@credits}{} % If not set in options
\newcommand{\liumaiex@regnumber}{}
\newcommand{\liumaiex@svkeywords}{}
\newcommand{\liumaiex@enkeywords}{}
\newcommand{\liumaiex@publishurl}{}
\newcommand{\liumaiex@metadataauthor}{\liumaiex@author}
\newcommand{\liumaiex@metadatatitle}{\liumaiex@title}
\newcommand{\liumaiex@metadatakeywords}{}
\newcommand{\liumaiex@metadatasubject}{}

% Labels, will be set according to language
\newcommand{\liu}{}
\newcommand{\mai}{}
\newcommand{\liumaiex@creditslabel}{}
\newcommand{\liumaiex@levellabel}{}
\newcommand{\liumaiex@supervisorlabel}{}
\newcommand{\liumaiex@multisupervisorlabel}{}
\newcommand{\liumaiex@examinerlabel}{}
\newcommand{\liumaiex@svkeywordslabel}{Nyckelord}
\newcommand{\liumaiex@enkeywordslabel}{Keywords}
\newcommand{\liumaiex@svurllabel}{URL för elektronisk version}
\newcommand{\liumaiex@enurllabel}{URL for electronic version}


%---------------------------------------
%	Language
%---------------------------------------

% liumaiex@langswe set in options
% Swedish and English loaded by Babel/Polyglossia in packages
\ifliumaiex@langswe
	\renewcommand{\liu}{Linköpings universitet}
	\renewcommand{\mai}{Matematiska institutionen, \liu}
	\renewcommand{\liumaiex@creditslabel}{Högskolepoäng}
	\renewcommand{\liumaiex@levellabel}{Nivå}
	\renewcommand{\liumaiex@supervisorlabel}{Handledare}
	\renewcommand{\liumaiex@multisupervisorlabel}{Handledare}
	\renewcommand{\liumaiex@examinerlabel}{Examinator}
	\renewcommand{\liumaiex@metadatakeywords}{\liumaiex@svkeywords}
	\ifliumaiex@xeorlua
		\setmainlanguage{swedish}
		\setotherlanguage{english}
	\else
		\AtBeginDocument{\selectlanguage{swedish}}
	\fi
\else
	\renewcommand{\liu}{Linköping University}
	\renewcommand{\mai}{Department of Mathematics, \liu}
	\renewcommand{\liumaiex@creditslabel}{Credits}
	\renewcommand{\liumaiex@levellabel}{Level}
	\renewcommand{\liumaiex@supervisorlabel}{Supervisor}
	\renewcommand{\liumaiex@multisupervisorlabel}{Supervisors}
	\renewcommand{\liumaiex@examinerlabel}{Examiner}
	\renewcommand{\liumaiex@metadatakeywords}{\liumaiex@enkeywords}
	\ifliumaiex@xeorlua
		\setmainlanguage{english}
		\setotherlanguage{swedish}
	\else
		\AtBeginDocument{\selectlanguage{english}}
	\fi
\fi


%---------------------------------------
%	Info commands
%---------------------------------------

\renewcommand*{\title}[1]{\renewcommand{\liumaiex@title}{#1}}
\renewcommand*{\author}[1]{\renewcommand{\liumaiex@author}{#1}}
\newcommand*{\shortauthor}[1]{\renewcommand{\liumaiex@shortauthor}{#1}}
\newcommand*{\publishmonth}[1]{\renewcommand{\liumaiex@publishmonth}{#1}}
\newcommand*{\publishyear}[1]{\renewcommand{\liumaiex@publishyear}{#1}}
\newcommand*{\city}[1]{\renewcommand{\liumaiex@city}{#1}}
\newcommand*{\department}[1]{\renewcommand{\liumaiex@department}{#1}}
\newcommand*{\level}[1]{\renewcommand{\liumaiex@level}{#1}}
\newcommand*{\credits}[1]{\renewcommand{\liumaiex@credits}{#1}}
\newcommand*{\regnumber}[1]{\renewcommand{\liumaiex@regnumber}{#1}}
\newcommand*{\svkeywords}[1]{\renewcommand{\liumaiex@svkeywords}{#1}}
\newcommand*{\enkeywords}[1]{\renewcommand{\liumaiex@enkeywords}{#1}}
\newcommand*{\publishurl}[1]{\renewcommand{\liumaiex@publishurl}{#1}}
\newcommand*{\pdfauthor}[1]{\renewcommand{\liumaiex@metadataauthor}{#1}}
\newcommand*{\pdftitle}[1]{\renewcommand{\liumaiex@metadatatitle}{#1}}
\newcommand*{\pdfkeywords}[1]{\renewcommand{\liumaiex@metadatakeywords}{#1}}
\newcommand*{\pdfsubject}[1]{\renewcommand{\liumaiex@metadatasubject}{#1}}

% The following code allows examiner to be entered with or without star
\newcommand{\examiner}
	{\@ifstar\liumaiex@examiner@star\liumaiex@examiner@nostar}
\newcommand*{\liumaiex@examiner@star}[1]
	{\liumaiex@examiner@nostar{#1}{\liumaiex@department}}
\newcommand*{\liumaiex@examiner@nostar}[2]{%
	\renewcommand*{\liumaiex@examiner}{{#1}{#2}}}

% The following code allows multiple supervisors to be entered
% (with or without star). They are stored as e.g. \liumaiex@supervisor1,
% and the information is stored as {<name>}{<department>}
\newcounter{liumaiex@numberofsupervisors}
\newcommand{\supervisor}
	{\@ifstar\liumaiex@supervisor@star\liumaiex@supervisor@nostar}
\newcommand*{\liumaiex@supervisor@star}[1]
	{\liumaiex@supervisor@nostar{#1}{\liumaiex@department}}
\newcommand*{\liumaiex@supervisor@nostar}[2]{%
	\stepcounter{liumaiex@numberofsupervisors}
	\expandafter\newcommand
		\csname liumaiex@supervisor\the\value{liumaiex@numberofsupervisors}%
			\endcsname
		{{#1}{#2}}}


%---------------------------------------
%	Header
%---------------------------------------

\fancypagestyle{normal}{
	\fancyhf{}
	\fancyhead[LO]{\nouppercase{\small\rightmark}}
	\fancyhead[RE]{\nouppercase{\small\leftmark}}
	\fancyhead[LE,RO]{\small\thepage}
	\renewcommand{\headrulewidth}{0.2pt}
	\renewcommand{\footrulewidth}{0pt}}
\fancypagestyle{plain}{
	\fancyhf{}
	\fancyfoot[RE,LO]{\small\liumaiex@shortauthor,\space\liumaiex@publishyear.}
	\fancyfoot[LE,RO]{\small\thepage}
	\renewcommand{\headrulewidth}{0pt}
	\renewcommand{\footrulewidth}{0.2pt}}
\fancypagestyle{copyright}{
	\fancyhf{}
	\fancyfoot[L]{\copyright\space\liumaiex@publishyear,\space\liumaiex@author}
	\renewcommand{\headrulewidth}{0pt}}
\pagestyle{normal}

% This redefines \cleardoublepage to remove headers between chapters
\def\cleardoublepage{%
	\clearpage
	% \if@twoside
	% 	\ifodd\c@page
	% 		% do nothing
	% 	\else
	% 		\hbox{}%
	% 		\thispagestyle{empty}%
	% 		\newpage
	% 		\if@twocolumn
	% 			\hbox{}%
	% 			\newpage
	% 		\fi
	% 	\fi
	% \fi
	}


%---------------------------------------
%	Title
%---------------------------------------

\newcommand*{\liumaiex@supervisorexaminerformat}[2]{\textbf{#1},\\#2}
\newcounter{liumaiex@currentsupervisor}
\newcommand{\liumaiex@endline}{\rule[-1em]{0pt}{1em}\\}
% \liumaiex@endline is a workaround to get spacing to work
% after minipages, replace with \rule[-1em]{1pt}{1em}\\
% to see how it works
\renewcommand{\maketitle}{%
	\cleardoublepage
	\thispagestyle{empty}%
	\begin{center}
		\null
		\vskip 4cm%
		{\Large\textbf{\liumaiex@title}\vskip 1em}%
		\liumaiex@department\vskip 1em%
		\textbf{\liumaiex@author}\vskip 1em%
		\liumaiex@regnumber
	\end{center}
	\vfill
	\begin{tabularx}{\textwidth}{@{}r@{\quad}X@{}}
		\liumaiex@creditslabel:&
			\textbf{\liumaiex@credits}\liumaiex@endline
		\liumaiex@levellabel:&
			\textbf{\liumaiex@level}\liumaiex@endline
		\ifnum\value{liumaiex@numberofsupervisors}>1\relax
			\liumaiex@multisupervisorlabel:
		\else
			\liumaiex@supervisorlabel:
		\fi&
			\begin{minipage}[t]{\linewidth}\raggedright
				\loop\ifnum\value{liumaiex@currentsupervisor}<
						\value{liumaiex@numberofsupervisors}%
					\stepcounter{liumaiex@currentsupervisor}%
					\expandafter\expandafter\expandafter
						\liumaiex@supervisorexaminerformat
							\csname liumaiex@supervisor%
								\the\value{liumaiex@currentsupervisor}%
								\endcsname\liumaiex@endline
					% the first and third \expandafter replaces
					% \csname...\endcsname with e.g.
					% \liumaiex@supervisor1, and the second \expandafter
					% replaces \liumaiex@supervisor1 with {...}{...},
					% so it becomes
					% \liumaiex@supervisorexaminerformat{...}{...}
				\repeat
			\end{minipage}\\
		\liumaiex@examinerlabel:&
			\begin{minipage}[t]{\linewidth}\raggedright
				\expandafter\liumaiex@supervisorexaminerformat\liumaiex@examiner
				\liumaiex@endline
			\end{minipage}\\
		\liumaiex@city:&
			\textbf{\liumaiex@publishmonth\space\liumaiex@publishyear}
	\end{tabularx}
	\liumaiex@placemetadata
	\cleardoublepage}


%---------------------------------------
%	After abstract
%---------------------------------------

\newcommand{\placesvkeywords}
	{\begin{description}
		\item[\liumaiex@svkeywordslabel:]\leavevmode\\
		\liumaiex@svkeywords
	\end{description}}
\newcommand{\placeenkeywords}
	{\begin{description}
		\item[\liumaiex@enkeywordslabel:]\leavevmode\\
		\liumaiex@enkeywords
	\end{description}}
\newcommand{\placesvurl}
	{\begin{description}
		\item[\liumaiex@svurllabel:]\leavevmode\\
		\url{\liumaiex@publishurl}
	\end{description}}
\newcommand{\placeenurl}
	{\begin{description}
		\item[\liumaiex@enurllabel:]\leavevmode\\
		\url{\liumaiex@publishurl}
	\end{description}}


%---------------------------------------
%	Logo
%---------------------------------------

\expandafter\def\csname LiUMAIeX \endcsname{L\kern -.2em\raise .58ex\hbox {i}\kern -.07em\lower .25ex\hbox {U}\kern -.1emM\kern -.12em{\sbox \z@ T\vbox to\ht \z@ {\hbox {\check@mathfonts \fontsize \sf@size \z@ \math@fontsfalse \selectfont A}\vss }}\kern -.11emI\kern -.1em\lower .5ex\hbox {E}\kern -.125emX}%
\expandafter\def\expandafter\LiUMAIeX\expandafter{\expandafter\protect\csname LiUMAIeX \endcsname}%


%---------------------------------------
%	PDF metadata
%---------------------------------------

\newcommand{\liumaiex@placemetadata}{%
	\ifdefined\hypersetup % i.e. if hyperref is loaded
		\hypersetup{%
			pdftitle=\liumaiex@metadatatitle,
			pdfauthor=\liumaiex@metadataauthor,
			pdfsubject=\liumaiex@metadatasubject,
			pdfkeywords=\liumaiex@metadatakeywords}
	\else
		\ifdefined\pdfinfo % works with pdftex and luatex
			\pdfinfo{
				/Title (\liumaiex@metadatatitle)
				/Author (\liumaiex@metadataauthor)
				/Keywords (\liumaiex@metadatakeywords)
				/Subject (\liumaiex@metadatasubject)}
		\fi
	\fi}


%---------------------------------------
%	Copyright page
%---------------------------------------

% This page is based on the template available from
% Linköping University Electronic Press,
% at the time of writing available at
% http://www.ep.liu.se/copyright/copyright_thesis.rtf
% (linked from http://www.ep.liu.se/exjobb/index.en.asp)
%
% The logo is (at the time of writing) available from
% http://www.liu.se/insidan/kommunikationsstod/grafiskprofil

% \AtEndDocument{%
% \clearpage\null\pagestyle{empty}\cleardoublepage
% \thispagestyle{copyright}
% \noindent
% {\large Linköping University Electronic Press}
% \hfill
% \ifliumaiex@langswe
% 	\smash{\raisebox{-9.7pt}{%
% 		\resizebox{5cm}{!}{\includegraphics{liulogo_sv}}}}%
% 	\hspace{-9.5pt}
% \else
% 	\smash{\raisebox{-10.2pt}{%
% 		\resizebox{5cm}{!}{\includegraphics{liulogo_en}}}}%
% 	\hspace{-10pt}
% \fi
% % The distances above are calibrated for the current logo;
% % if the logo changes, the distances probably should too.
% \par\vfil
% \bigskip
% \ifliumaiex@langswe
% 	\liumaiex@swedishcopyrighttext
% 	\liumaiex@englishcopyrighttext
% \else
% 	\liumaiex@englishcopyrighttext
% 	\liumaiex@swedishcopyrighttext
% \fi}

\newcommand{\liumaiex@swedishcopyrighttext}{%
\begin{otherlanguage*}{swedish}
\subsection*{Upphovsrätt}
Detta dokument hålls tillgängligt på Internet – eller dess framtida ersättare – från publiceringsdatum under förutsättning att inga extraordinära omständigheter uppstår.

Tillgång till dokumentet innebär tillstånd för var och en att läsa, ladda ner, skriva ut enstaka kopior för enskilt bruk och att använda det oförändrat för ickekommersiell forskning och för undervisning. Överföring av upphovsrätten vid en senare tidpunkt kan inte upphäva detta tillstånd. All annan användning av dokumentet kräver upphovsmannens medgivande. För att garantera äktheten, säkerheten och tillgängligheten finns lösningar av teknisk och administrativ art.

Upphovsmannens ideella rätt innefattar rätt att bli nämnd som upphovsman i den omfattning som god sed kräver vid användning av dokumentet på ovan beskrivna sätt samt skydd mot att dokumentet ändras eller presenteras i sådan form eller i sådant sammanhang som är kränkande för upphovsmannens litterära eller konstnärliga anseende eller egenart.

För ytterligare information om Linköping University Electronic Press se förlagets hemsida \url{http://www.ep.liu.se/}.
\end{otherlanguage*}}

\newcommand{\liumaiex@englishcopyrighttext}{%
\begin{otherlanguage*}{english}
\subsection*{Copyright}
The publishers will keep this document online on the Internet – or its possible replacement – from the date of publication barring exceptional circumstances.

The online availability of the document implies permanent permission for anyone to read, to download, or to print out single copies for his/her own use and to use it unchanged for non-commercial research and educational purpose. Subsequent transfers of copyright cannot revoke this permission. All other uses of the document are conditional upon the consent of the copyright owner. The publisher has taken technical and administrative measures to assure authenticity, security and accessibility.

According to intellectual property law the author has the right to be mentioned when his/her work is accessed as described above and to be protected against infringement.

For additional information about the Linköping University Electronic Press and its procedures for publication and for assurance of document integrity, please refer to its www home page: \url{http://www.ep.liu.se/}.
\end{otherlanguage*}}


\endinput
