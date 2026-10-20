UML-диаграмма классов
=====================

Диаграмма классов показывает отношения наследования между основными компонентами программы. Стрелка с незакрашенным наконечником обозначает обобщение: стрелка направлена к базовому классу.

.. graphviz::
   :caption: UML-диаграмма классов проекта

   digraph UML {
       graph [rankdir=TB, bgcolor="transparent", pad="0.25", nodesep="0.42", ranksep="0.72", splines=ortho];
       node [shape=plain, fontname="DejaVu Sans", fontsize=10];
       edge [color="#477653", penwidth=1.5, arrowsize=1.0, arrowhead=empty];

       Reader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#B6E3C8"><I>«abstract»</I><BR/><B>Reader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Общий интерфейс чтения</TD></TR>
         </TABLE>>];

       SequenceReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#CBEBD5"><B>SequenceReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Работа с последовательностями</TD></TR>
         </TABLE>>];

       GenomicDataReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#B6E3C8"><I>«abstract»</I><BR/><B>GenomicDataReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Обработка геномных записей</TD></TR>
         </TABLE>>];

       FastaReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#DDF2E3"><B>FastaReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Формат FASTA</TD></TR>
         </TABLE>>];

       FastqReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#DDF2E3"><B>FastqReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Формат FASTQ</TD></TR>
         </TABLE>>];

       SamReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#DDF2E3"><B>SamReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Формат SAM</TD></TR>
         </TABLE>>];

       VcfReader [label=<
         <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="8" COLOR="#548767">
           <TR><TD BGCOLOR="#DDF2E3"><B>VcfReader</B></TD></TR>
           <TR><TD BGCOLOR="#F6FCF7">Формат VCF</TD></TR>
         </TABLE>>];

       SequenceReader -> Reader;
       FastaReader -> SequenceReader;
       FastqReader -> SequenceReader;
       SamReader -> GenomicDataReader;
       VcfReader -> GenomicDataReader;

       {rank=same; SequenceReader; GenomicDataReader}
       {rank=same; FastaReader; FastqReader; SamReader; VcfReader}
   }

Классы сгруппированы по назначению: обработка последовательностей и обработка геномных записей.
