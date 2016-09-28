package ara.dr

import java.io.File
import java.io.PrintWriter
import atk.compbio.gff._
import atk.util.Tool
import ara.dr.DRVcfLine
import ara.dr.DRsnp._
import java.nio.charset.CodingErrorAction
import scala.io.Codec
import ara.CodonConfig
import ara.MacawSNPtyper
import scala.io.Codec.decoder2codec
import ara.dr.DRVcfLine

/**
 *  Read VCF-files of samples mapped against a minimized reference genome,
 *  consisting of only drug resistance regions.
 *
 */

object DrugResistanceSnps extends CodonConfig with Tool {

  case class Config(val vcf: File = null, val gff: File = null, val fasta: File = null, val output: String = null, val all: Boolean = false)

  def complement(seq: String) = {
    seq.map { c =>
      c match {
        case 'A' => 'T'
        case 'T' => 'A'
        case 'C' => 'G'
        case 'G' => 'C'
        case _ => 'N'
      }
    }.mkString
  }

  val decoder = Codec.UTF8.decoder.onMalformedInput(CodingErrorAction.IGNORE)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar dr-snps") {
      opt[File]("vcf") required () action { (x, c) => c.copy(vcf = x) } text ("VCF-file")
      opt[File]("gff") required () action { (x, c) => c.copy(gff = x) } text ("GFF-file")
      opt[File]("fasta") required () action { (x, c) => c.copy(fasta = x) } text ("Fasta-file")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x + ".dr-snps.txt") } text ("Output name.")
      opt[Unit]("all") action { (x, c) => c.copy(all = true) }  text ("Output all mutations as comments for debug purposes")
    }

    /* Load Coll2015 library of drug resistance SNPs */
    val drList = scala.io.Source.fromInputStream(MacawSNPtyper.getClass().getResourceAsStream("/Coll2015DrugResistances.txt"))(decoder).getLines().filterNot(_.startsWith("#")).filter(_.isSNP).map { line =>
      line match {
        case DRsnp(d, l, lt, cp, r, gp, a, cn, aac) => new DRsnp(d, l, lt, cp, r, gp, a, cn, aac)
      }
    }.toList
    val locusTag = drList.map(d => ((if (d.locus.endsWith("-promoter")) d.locus.dropRight(9) else d.locus) -> d.locusTag)).toMap
    val associatedDrug = drList.map(d => ((if (d.locus.endsWith("-promoter")) d.locus.dropRight(9) else d.locus) -> d.drug)).groupBy(_._1).mapValues(_.map(_._2).distinct.mkString(","))

    /* Find gene region of SNP position */
    def getLocus(cp: Int, lociList: List[(String, GFFLine)], ref: String, alt: String): List[Map[String, Any]] = {
      val loci = lociList.sortBy(_._2.start)

      val l = (0 until loci.size).map { idx =>
        val locus = loci(idx)._1
        val drug = associatedDrug(locus)
        val strand = loci(idx)._2.line.split("\t")(6)
        val start = loci(idx)._2.start
        val end = loci(idx)._2.end
        if (strand.equals("+")) {
          if (idx == 0) {
            if (cp < start) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (cp - start), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (cp - start + 1), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else if (loci.size == 1 && cp > end) Map("locus" -> (locus + "-tail"), "drug" -> drug, "gene-coordinate" -> (end - cp), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else (null)
          } else if (idx == loci.size - 1) {
            if (cp < start && cp > loci(idx - 1)._2.end) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (cp - start), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (cp - start + 1), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else if (cp > end) Map("locus" -> (locus + "-tail"), "drug" -> drug, "gene-coordinate" -> (cp - start + 1), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else (null)
          } else {
            if (cp < start && cp > loci(idx - 1)._2.end) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (cp - start), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (cp - start + 1), "nucleotide-change" -> (ref + "/" + alt), "strand" -> "+")
            else (null)
          }
        } else { //negative strand
          if (idx == loci.size - 1) {
            if (cp > end) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (end - cp), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (end - cp + 1), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else if (cp < start && loci.size == 1) Map("locus" -> (locus + "-tail"), "drug" -> drug, "gene-coordinate" -> (end - cp + 1), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else (null)
          } else if (idx == 0) {
            if (cp < start) Map("locus" -> (locus + "-tail"), "drug" -> drug, "gene-coordinate" -> (end - cp + 1), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (end - cp + 1), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else if (cp > end && cp < loci(idx + 1)._2.start) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (end - cp), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else (null)
          } else {
            if (cp > end && cp < loci(idx + 1)._2.start) Map("locus" -> (locus + "-promoter"), "drug" -> drug, "gene-coordinate" -> (end - cp), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else if (cp >= start && cp <= end) Map("locus" -> locus, "drug" -> drug, "gene-coordinate" -> (end - cp + 1), "nucleotide-change" -> (complement(ref) + "/" + complement(alt)), "strand" -> "-")
            else (null)
          }
        }
      }.toList
      l.filterNot(_ == null)
    }

    parser.parse(args, Config()) map { config =>

      /* Prepare GFF map */
      val gff = GFFFile(config.gff).filterNot(_.kind.equals("CDS")).filterNot(_.kind.equals("source")).filterNot(_.line.split("\t")(8).startsWith("note="))
      val gffGenes = gff.map(g => (g.attributes("locus_tag").split(""""""")(1) -> g)).toMap

      /* Load full reference genome fasta file */
      val ref = tLines(config.fasta).filterNot(_.startsWith(">")).mkString

      /* Read VCF file and filter SNPs */
      val snps = tLines(config.vcf).map(DRVcfLine(_)).filter(_.isSNP())

      val snpInfo = snps.map { snp =>
        val loci = snp.loci.map { locus => (locus, gffGenes(locusTag(locus))) }.toList
        val locus = getLocus(snp.chrPos, loci, snp.ref, snp.alt)

        val snpInfoPerLocus = locus.map { l =>

          val (codonNumber, codonChange, aminoAcidChange, knownInfo) = { //codonChange on forward strand
            val locusName = l("locus").asInstanceOf[String]
            if (locusName.endsWith("-tail") || locusName.endsWith("-promoter") || locusName.equals("rrl") || locusName.equals("rrs")) {
              val knownChrPos = drList.filter(_.cp.equals(snp.chrPos))
              val knownSnp = if (knownChrPos.isEmpty) "Unknown mutation"
              else {
                val knownNchange = knownChrPos.filter(s => s.r.equals(snp.ref) && s.a.equals(snp.alt))
                if (knownNchange.isEmpty) "Known nucleotide coordinate"
                else "Known mutation"
              }
              ("-", "-", "-", knownSnp)
            } else {
              val gc = ((l("gene-coordinate")).asInstanceOf[Int])
              val cn = (gc + 2) / 3
              val codonPosition = (gc + 2) % 3 + 1
              var cc = { // codon change on forward strand
                if (codonPosition == 1) ref.substring(snp.chrPos - 1, snp.chrPos + 2) + "/" + snp.alt + ref.substring(snp.chrPos, snp.chrPos + 2)
                else if (codonPosition == 2) ref.substring(snp.chrPos - 2, snp.chrPos + 1) + "/" + ref.charAt(snp.chrPos - 2) + snp.alt + ref.charAt(snp.chrPos)
                else ref.substring(snp.chrPos - 3, snp.chrPos) + "/" + ref.substring(snp.chrPos - 3, snp.chrPos - 1) + snp.alt
              }
              if (gffGenes(locusTag(l("locus").asInstanceOf[String])).line.split("\t")(6).equals("-")) cc = cc.split("/").map(c => complement(c.reverse)).mkString("/")
              val ac = cc.split("/").map(codonMap(_)).mkString("/")
              val knownCodon = {
                val knownCodonNumber = drList.filter(s => s.locus.equals(l("locus")) && s.gp == gc)
                if (knownCodonNumber.isEmpty) "Unknown mutation"
                else {
                  val knownAac = knownCodonNumber.filter(_.aaChange.equals(ac))
                  if (knownAac.isEmpty) "Known amino acid coordinate."
                  else "Known mutation"
                }
              }
              (cn, cc, ac, knownCodon)
            }
          }

          val nCount = if (l("strand") == "+") {
            snp.bc.drop(3)
          } else {
            val acgt = snp.bc.drop(3).split(",")
            acgt.reverse.mkString(",")
          }

          (Map("drug" -> l("drug"), "locus" -> l("locus"), "chrPos" -> snp.chrPos, "gene-coordinate" -> l("gene-coordinate"), "nucleotide-change" -> l("nucleotide-change"), "codon-number" -> codonNumber, "codon-change" -> codonChange, "aa-change" -> aminoAcidChange, "mut-info" -> knownInfo, "filter" -> snp.filter, "count" -> nCount))
        }
        snpInfoPerLocus
      }.flatten

      val knownMutations = snpInfo.filter(m => m("mut-info") == "Known mutation")

      val pw = new PrintWriter(config.output)
      pw.println("# SNPs: " + snps.size)
      pw.println("# Known SNPs: " + knownMutations.size)
      pw.println("# Ambiguous known SNPs: " + knownMutations.filter(m => m("filter") == "Amb").size)
      pw.println("#")
      pw.println("# Known mutations in drug resistance regions")
      pw.println("# Drug\tLocus\tChromosome coordinate\tGene coordinate\tNucleotide change\tCodon number\tCodon change\tAmino acid change\tMutation info\tFilter\tCount")
      knownMutations.foreach { l =>
        pw.println(l("drug") + "\t" + l("locus") + "\t" + l("chrPos") + "\t" + l("gene-coordinate") + "\t" + l("nucleotide-change") + "\t" + l("codon-number") + "\t" + l("codon-change") + "\t" + l("aa-change") + "\t" + l("mut-info") + "\t" + l("filter") + "\t" + l("count"))
      }
      pw.println("## All SNPS")
      pw.println("##----------")
      for (snp <- snps) {

        pw.println("# " + snp.line)
      }
      pw.close

    }

  }
}