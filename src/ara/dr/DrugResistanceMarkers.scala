package ara.dr

import scala.io.Source
import java.io.File
import java.io.PrintWriter
import ara.Gene._
import java.util.Calendar
import ara.CodonConfig
import ara.Gene

/**
 * Reads only SNPs from drug resistance list, and produces markers
 */

object DrugResistanceMarkers extends CodonConfig{

  case class Config(
    val drList: File = null,
    val fasta: File = null,
    val output: String = null,
    val gff: File = null,
    val supportingRefs: Int = 5
    )

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar dr-markers") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(drList = x) } text ("List of TB drug resistances.")
      opt[File]("ref") required () action { (x, c) => c.copy(fasta = x) } text ("Reference fasta file.")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Name of output file.")
      opt[File]("gff") required () action { (x, c) => c.copy(gff = x) } text ("gff-file corresponding to reference fasta-file.")
      opt[Int]("supRef") action { (x, c) => c.copy(supportingRefs = x) } text ("Threshold to include SNPs based on number of supporting references. Default = 5")
      }

    parser.parse(args, Config()) map { config =>

      val ref = Source.fromFile(config.fasta).getLines.filterNot(_.startsWith(">")).mkString
      val supRef = config.supportingRefs
      
      class DRsnp(val drug: String, val locus: String, val locusTag: String, val r: String, val p: Int, val a: String, val cause: String) extends Ordered[DRsnp] {
        override def toString(): String = ">" + r + p + a + "_" + locus + "_" + cause + "_" + drug
        def complementToString(drugs: String): String = ">" + r + p + r + "_susceptibility_" + drugs
        def compare(that: DRsnp): Int = this.p compare that.p
        val marker = ref.substring(p - 11, p - 1) + a + ref.substring(p, p + 10)
      }

      object DRsnp {
        def unapply(s: String): Option[(String, String, String, String, Int, String)] = {
          val sArr = s.mkString.split("\t")
          val chrPos = sArr(2)
          if (!chrPos.contains("/") && !chrPos.equals("-")) {
            val drug = sArr(0)
            val locus = if (sArr(1).contains("_")) sArr(1).split("_").mkString("-") else sArr(1)
            val locusTag = sArr(8)
            val totalSources = sArr(9).split(";").size
            val nChange = sArr(4)
            val ncArr = nChange.split("/")
            val r = ncArr(0)
            val a = ncArr(1)
            val nucleotides = Array[String]("A", "C", "T", "G")
            if (nucleotides.contains(r) && nucleotides.contains(a) && totalSources >= supRef)
              Some((drug, locus, locusTag, r, chrPos.toInt, a))
            else None
          } else None
        }
      }

      class Mutation(s: String) {
        def isSNP: Boolean = s match {
          case DRsnp(d, l, lt, r, p, a) => true
          case _ => false
        }
      }

      implicit def seqtoBool(s: String) = new Mutation(s)

      class compSeq(val seq: String) {
        def complement = {
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
      }

      implicit def seqToSeq(s: String) = new compSeq(s)



      val drList = Source.fromFile(config.drList).getLines.filterNot(_.startsWith("#")).filter(_.isSNP).map { line =>
        line match {
          case DRsnp(d, l, lt, r, p, a) => new DRsnp(d, l, lt, r, p, a, "resistance")
        }
      }.toList

      /** SNPs grouped by position and mutation.*/
      val snpsPerPosition = drList.groupBy(_.p).map(snp => snp match {
        case (p, snps) => (p, snps.groupBy(_.a))
      }).toList.sortBy(_._1)

      println(snpsPerPosition.flatMap(_._2.flatMap(_._2)).size + " SNPs")

      /** Unique SNPs per position. Duplicate mutations indicating resistance to other drugs are merged. */
      val uniqueSnps = snpsPerPosition.map(pos => pos match {
        case (p, snpMap) => (p, snpMap.map(mutation => mutation match {
          case (alt, snps) => snps.size match {
            case x if x > 1 =>
              val drugs = snps.map(_.drug).mkString("_"); (alt, List(new DRsnp(drugs, snps.head.locus, snps.head.locusTag, snps.head.r, p, alt, "resistance")))
            case 1 => (alt, snps)
          }
        }).flatMap(_._2).toList)
      })

      println(uniqueSnps.flatMap(_._2).size + " unique SNPs")
      println(uniqueSnps.size + " unique SNP coordinations")

      /**
       * val gff = Source.fromFile(config.gff).getLines.filterNot(line => line.split("\t")(2) == "CDS" && line.split("\t")(2) == "source").map(line => line match {
       * case Gene(tag, start, end, dir) => new Gene(tag, start, end, dir)
       * }).toList
       */

      val gff = Source.fromFile(config.gff).getLines.filterNot { line =>
        val feature = line.split("\t")(2)
        feature == "CDS" || feature == "source"
      }.map {
        _ match {
          case Gene(tag, feature, start, end, dir) => (tag -> new Gene(tag, feature, start, end, dir))
        }
      }.toMap.filter {
        _ match {
          case (tag, gene) => uniqueSnps.map(_._2.head.locusTag).contains(tag)
        }
      }
      //gff.foreach(println)

      def mutationEffect(unknownBases: List[DRsnp], gene: Gene, pos: Int): List[(DRsnp, String)] = {
        //val pos = unknownBases.head.p
        val codonCoordination = if (gene.dir == "+") (pos - gene.start + 1) % 3 else (gene.end - pos + 1) % 3
        val codonChange = codonCoordination match {
          case 0 => if (gene.dir == "+") unknownBases.map(snp => (snp, ref.substring(pos - 3, pos - 1) + snp.a)) else unknownBases.map(snp => (snp, (snp.a + ref.substring(pos, pos + 2)).reverse.complement))
          case 1 => if (gene.dir == "+") unknownBases.map(snp => (snp, snp.a + ref.substring(pos, pos + 2))) else unknownBases.map(snp => (snp, (ref.substring(pos - 3, pos - 1) + snp.a).reverse.complement))
          case 2 => {
            val codons = unknownBases.map(snp => (snp, ref.charAt(pos - 2) + snp.a + ref.charAt(pos)))
            if (gene.dir == "+") codons
            else codons.map(ac => (ac._1, ac._2.reverse.complement))
          }
        }
        val refCodon = codonCoordination match {
          case 0 => if (gene.dir == "+") ref.substring(pos - 3, pos) else ref.substring(pos - 1, pos + 2).reverse.complement
          case 1 => if (gene.dir == "+") ref.substring(pos - 1, pos + 2) else ref.substring(pos - 3, pos).reverse.complement
          case 2 => {
            val codon = ref.substring(pos - 2, pos + 1)
            if (gene.dir == "+") codon
            else codon.reverse.complement
          }
        }
        val refAA = codonMap(refCodon)
        val AAChange = codonChange.map { ac =>
          if (codonMap(ac._2) == refAA) (ac._1, ac._2, codonMap(ac._2), "synonymous")
          else (ac._1, ac._2, codonMap(ac._2), "nonsynonymous")
        }
        // println("Refcodon:\t" + refCodon + "\tRefAA:\t" + refAA + "\tAAChange:\t" + AAChange)
        AAChange.map(snp => (snp._1, snp._4))
      }

      /** Convert unique SNPs to markers, and add susceptible marker per genome position. */
      val markers = uniqueSnps.flatMap(_ match {
        case (pos, snps) =>
          //println(snps.size + " snps: " + snps) 
          val alts = snps.map(_.a)
          val r = snps.head.r
          val locus = snps.head.locus
          val locusTag = snps.head.locusTag
          val unknownBases = List("A", "C", "G", "T").filter(!_.equals(r)).filter(!alts.contains(_)) //.map(new DRsnp("unknown", locus, locusTag, r, pos, _, "unknown"))
          val gene = gff(locusTag)

          val unknownSnps = unknownBases.map(alt => new DRsnp("unknown", locus, locusTag, r, pos, alt, "unknown"))

          val unknownSnps2 = if (gene.feature.endsWith("RNA") || locus.endsWith("promoter"))
            unknownSnps.map { snp => (">" + snp.r + snp.p + snp.a + "_" + snp.locus + "_" + "-" + "_" + snp.cause + "_" + snp.drug, snp.marker) }
          else {
            mutationEffect(unknownSnps, gene, pos).map {
              _ match {
                case (snp, ns) => (">" + snp.r + snp.p + snp.a + "_" + snp.locus + "_" + ns + "_" + snp.cause + "_" + snp.drug, snp.marker)
              }
            }
          }

          val knownSnps = snps.map { snp =>
            (">" + snp.r + snp.p + snp.a + "_" + snp.locus + "_" + "nonsynonymous" + "_" + snp.cause + "_" + snp.drug, snp.marker)
          }
          val mList = unknownSnps2 ++ knownSnps
          val drugs = snps.flatMap(d => if (d.drug.contains("_")) d.drug.split("_").toList else List(d.drug)).distinct.mkString("_")
          mList :+ (">" + r + pos + r + "_" + locus + "_" + "-" + "_susceptibility_" + drugs, ref.substring(pos - 11, pos + 10))
      })

      //markers.foreach(println)
      println(markers.size + " markers")

      /**
       * Filter out non-unique markers
       */
      val uniqueMarkers = markers.map(_ match {
        case (id, marker) => (id, marker, markers.count(_ match {
          case (id2, marker2) => marker == marker2
        }))
      }).filter(_._3 == 1).map(m => (m._1, m._2))

      println(uniqueMarkers.size + " unique markers")

      /** Print unique markers to output file. */
      val pw = new PrintWriter(new File(config.output))
      pw.println("# Drug resistance markers")
      pw.println("# SNPs with at least " + supRef + " supporting sources.")
      pw.println("# Compiled: " + Calendar.getInstance.getTime)
      markers.foreach(_ match {
        case (id, marker) => pw.println(id + "\n" + marker)
      })
      pw.close

    }

  }
}