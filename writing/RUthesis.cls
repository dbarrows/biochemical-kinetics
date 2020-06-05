%%
%% This is file `RUthesis.cls'.
%%
% Copyright 2010 by Catherine Beauchemin.
% Based on earlier work by Kipp Cannon which contained the
% following notice:
% "Copyright 2003 by Kipp Cannon."
% Based on earlier work by R. H. A. David Shaw which contained the
% following notice:
% "Copyright 2000 by R. H. A. David Shaw.
% This class file may be freely copied, distributed and modified with the
% sole proviso that this copyright notice be included.  Neither this style
% file nor any modification of it may be sold for profit."
%%
%% This file is designed for use with the LaTeX2e document formatting
%% system.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{RUthesis}[2010/01/30 v1.0 Ryerson University thesis class]
\RequirePackage{graphics}	% for including the RU crest
\RequirePackage{setspace}
\RequirePackage{fancyhdr}


%
% Process our class options.
%

\newif\if@RUgrad
\DeclareOption{grad}{\@RUgradtrue}
\DeclareOption{ugrad}{\@RUgradfalse}

\newif\if@RUlot
\DeclareOption{lot}{\@RUlottrue}
\DeclareOption{nolot}{\@RUlotfalse}

\newif\if@RUlof
\DeclareOption{lof}{\@RUloftrue}
\DeclareOption{nolof}{\@RUloffalse}

\newif\if@RUopenright
\DeclareOption{openright}{\@RUopenrighttrue}
\DeclareOption{openany}{\@RUopenrightfalse}

\newcommand{\@RUptsize}{}
\DeclareOption{10pt}{\renewcommand{\@RUptsize}{10pt}}
\DeclareOption{11pt}{\renewcommand{\@RUptsize}{11pt}}
\DeclareOption{12pt}{\renewcommand{\@RUptsize}{12pt}}

\newif\if@RUtwoside
\DeclareOption{oneside}{\@RUtwosidefalse}
\DeclareOption{twoside}{\@RUtwosidetrue}

\newcommand{\@RUspacing}{}
\newcommand{\@RUspacingenv}{}
\DeclareOption{singlespace}{\renewcommand{\@RUspacing}{\singlespacing}\renewcommand{\@RUspacingenv}{singlespace}}
\DeclareOption{onehalfspace}{\renewcommand{\@RUspacing}{\onehalfspacing}\renewcommand{\@RUspacingenv}{onehalfspace}}
\DeclareOption{doublespace}{\renewcommand{\@RUspacing}{\doublespacing}\renewcommand{\@RUspacingenv}{doublespace}}


\ExecuteOptions{nocrest,grad,lot,lof,openright,11pt,oneside,onehalfspace}
\ProcessOptions*


%
% Load the class upon which we are based.
%
\LoadClass[\if@RUopenright{openright}\else{openany}\fi,\@RUptsize,\if@RUtwoside{twoside}\else{oneside}\fi]{book}


%
% Fix the page numbering so it complies with SGS regulations.
%
\fancyhead{} % Clear all header fields
\if@RUtwoside{
	\fancyhead[RO]{\slshape \rightmark}
	\fancyhead[LE]{\slshape \leftmark}
}
\else{
	\fancyhead[R]{\slshape \leftmark}
    \fancyhead[L]{\slshape B.SC. THESIS - DEXTER BARROWS}
}\fi
\fancyfoot{} % Clear all footer fields
\fancyfoot[C]{\thepage}

%
% Set the names of some sections.  You can use \renewcommand{}{} in your
% preamble to redefine these, if needed (eg. if your thesis is in French).
%
\newcommand{\abstractname}{Abstract}
\newcommand{\acknowledgementsname}{Acknowledgements}
\renewcommand{\contentsname}{Table of Contents}
\newcommand{\dedicationpretext}{To}
\newcommand{\prefacename}{Preface}


%
% Set and adjust various lengths.
%
% These commands set the page margins just slightly larger than the FGSR
% regulation minimum.  The top, left, and right margins are approximately
% 3/32" larger than they need to be which should be enough to accomodate
% any potential page misalignment in the printer.  The bottom margin is a
% 1/4" larger than it needs to be but this makes the page look nicer than
% if things are crowded down as low as they can go.
%
\setlength{\oddsidemargin}{0.6in}
\setlength{\evensidemargin}{0.1in}
\setlength{\textwidth}{5.8in}
\setlength{\topmargin}{0in}
\setlength{\textheight}{8.25in}
\newlength{\RUsignatureheight}
\setlength{\RUsignatureheight}{36pt}


%
% This command is used to generate a place for a signature.
%
\newcommand{\RUsignature}[1]{
   \parbox{0.5\linewidth}{\vspace{\RUsignatureheight}\hrulefill\\#1}
}
\newcommand{\RUstudentsignature}{\hfill\RUsignature{\@author}\par}


%
% These commands generate the RUthesis class logo and version number.
%
\newcommand*{\RUthesis}{\texttt{RUthesis}}
\newcommand*{\RUthesisversion}{1.0}


%
% The required preamble commands are defined here.
%
\newcommand{\RUpastdegrees}[1]{%
   \newcommand{\theRUpastdegrees}{#1}%
}

\newcommand{\RUdegree}[1]{%
   \newcommand{\theRUdegree}{#1}%
}

\newcommand{\RUfield}[1]{%
   \newcommand{\theRUfield}{#1}%
}

\newcommand{\RUdepartment}[1]{%
   \newcommand{\theRUdept}{#1}%
}

\newcommand{\RUconvocation}[1]{%
   \newcommand{\theRUconv}{#1}%
}

\newcommand{\RUsubmitdate}[1]{%
   \newcommand{\theRUsubmitdate}{#1}%
}

\newcommand{\RUaddress}[1]{%
   \newcommand{\theRUaddress}{#1}%
}

\newcommand{\RUabstract}[1]{%
   \newcommand{\theRUabstract}{#1}%
}

\newcommand{\RUyear}[1]{%
   \newcommand{\theRUyear}{#1}%
}


%
% The optional preamble commands are defined here.
%
\newcommand{\RUsupervisor}[1]{%
   \newcommand{\theRUsupervisor}{#1}%
}

\newif\if@RUacknowledge
\@RUacknowledgefalse
\newcommand{\RUacknowledgement}[1]{%
   \@RUacknowledgetrue
   \newcommand{\theRUacknowledgement}{#1}%
}

\newif\if@RUdedication
\@RUdedicationfalse
\newcommand{\RUdedication}[1]{%
   \@RUdedicationtrue
   \newcommand{\theRUdedication}{#1}%
}

\newif\if@RUpreface
\@RUprefacefalse
\newcommand{\RUpreface}[1]{%
   \@RUprefacetrue
   \newcommand{\theRUpreface}{#1}%
}


%
% The prefatory pages are defined here.
%


%
% Frontispiece
%
\newcommand{\makeRUquote}[2]{
   \begin{quote}#1\end{quote}\par
   \begin{flushright}#2\makebox[1.0 in]{}\end{flushright}\par
}

\newcommand{\RUquotes}{}
\newif\if@RUfrontis
\@RUfrontisfalse
\newcommand{\RUaddquote}[2]{
   \@RUfrontistrue
   \toks0 = \expandafter{\RUquotes}
   \edef\RUquotes{\the\toks0 \noexpand\makeRUquote{#1}{#2}}
}

\newcommand{\makeRUfrontispiece}{
   \cleardoublepage
   \vspace*{\stretch{1}}
   \begin{\@RUspacingenv}
      \RUquotes
   \end{\@RUspacingenv}
   \vspace*{\stretch{2}}
   \@twosidetrue\clearpage\@twosidefalse
}


%
% Title Page
%
\newcommand{\makeRUtitlepage}{
	\cleardoublepage
	\begin{spacing}{2}
	\begin{center}
		\vspace{\stretch{3}}
			{\Large\MakeUppercase{\@title}}\par
		\vspace{\stretch{6}}
		by\par
		{\@author}\par
		\if@RUgrad{\theRUpastdegrees}\par
		\else{\theRUsubmitdate}\par\fi
		\vspace{\stretch{5}}
		{A thesis\\
		presented to \if@RUgrad{Ryerson University}
		\else{the Department of \theRUdept}\fi\\
		in partial fulfillment of the\\
		requirements for the degree of\\
		\theRUdegree \\
		in the Program of\\
		\theRUfield\\
		\if@RUgrad\else{at Ryerson University \\}\fi}\par
		\vspace{\stretch{20}}
		\if@RUgrad
		{Toronto, Ontario, Canada, \theRUyear}\par
		{\copyright\ \@author}
		\else
		{Supervisor: \theRUsupervisor}
		\fi
	\end{center}
	\end{spacing}
}


%
% Define the generic "prefpage" environment.
%
\newenvironment{RUprefpage}[1]{
   \cleardoublepage
   \centerline{\bfseries\large #1}\par
   \vspace{4ex}
}{
}

%
% Author Declaration Page
%
\newcommand{\makeRUauthordeclaration}{
	\begin{RUprefpage}{Author Declaration}
   		\vspace{2em}
		\noindent
		I hereby declare that I am the sole author of this thesis.\par
		\vspace{2em}
		\noindent
		I authorize Ryerson University to lend this thesis to other
		institutions or individuals for the purpose of scholarly
		research.\par
		\vspace{4em}
		\RUstudentsignature
		\vspace{10em}
		\noindent
		I further authorize Ryerson University to reproduce this thesis
		by photocopying or by other means, in total or in part, at the
		request of other institutions or individuals for the purpose
		of scholarly research.\par
		\vspace{4em}
		\RUstudentsignature
	\end{RUprefpage}
}


%
% Dedication
%
\newcommand{\makeRUdedicationpage}{
   \cleardoublepage
   \vspace{\stretch{1}}
   {\bfseries\itshape\Large
   \begin{center}
      \dedicationpretext\par
      \vspace{2ex}
      \theRUdedication
   \end{center}}
   \vspace*{\stretch{2.5}}
}


%
% Abstract
%
\newcommand{\makeRUabstractpage}{
	\begin{RUprefpage}{\abstractname}
		\if@RUgrad
		\begin{center}
			{\Large \@title}\par
			\vspace*{\stretch{0.5}}
			{\large \@author}\par
			\vspace{.5em}
			{\large \theRUdegree, \theRUfield}\par
			\vspace{.5em}
			{\large Ryerson University, \theRUyear}\par
			\vspace*{\stretch{1}}
		\end{center}
		\fi
	\end{RUprefpage}
	\doublespacing\theRUabstract
   \vspace{\stretch{5}}
}


%
% Preface
%
\newcommand{\makeRUpreface}{
   \begin{RUprefpage}{\prefacename}
      \begin{\@RUspacingenv}
         \theRUpreface
      \end{\@RUspacingenv}
   \end{RUprefpage}
}


%
% Acknowledgments
%
\newcommand{\makeRUacknowledgements}{
   \begin{RUprefpage}{\acknowledgementsname}
      \begin{\@RUspacingenv}
         \doublespacing\theRUacknowledgement
      \end{\@RUspacingenv}
   \end{RUprefpage}
}


%
% List of Tables
%
\newcommand{\makeRUlistoftables}{
   \cleardoublepage
   \begin{singlespace}
      \listoftables
   \end{singlespace}
   \cleardoublepage
}


%
% List of Figures
%
\newcommand{\makeRUlistoffigures}{
   \cleardoublepage
   \begin{singlespace}
      \listoffigures
   \end{singlespace}
   \cleardoublepage
}



%
% The maketitle command is redefined so as to generate the correct sequence
% of prefaratory pages.
%
\renewcommand{\maketitle}{
	\@twosidefalse
	\begin{titlepage}
		\if@RUfrontis\makeRUfrontispiece\fi
		\makeRUtitlepage
	\end{titlepage}
	\frontmatter
	\pagenumbering{roman}
	\pagestyle{plain}
	\setcounter{page}{2}
	\if@RUgrad\makeRUauthordeclaration\fi
	\makeRUabstractpage
	\if@RUpreface\makeRUpreface\fi
	\if@RUacknowledge\makeRUacknowledgements\fi
	\if@RUdedication\makeRUdedicationpage\fi
	\global\let\maketitle\relax
	\global\let\title\relax
	\global\let\author\relax
	\cleardoublepage
	\tableofcontents
	\if@RUlot\makeRUlistoftables\fi
	\if@RUlof\makeRUlistoffigures\fi
	\mainmatter
	\if@RUtwoside\@twosidetrue\fi
	\@RUspacing
	\pagestyle{fancy}
}
