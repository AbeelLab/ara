//package macaw2

import scala.io.Source
import java.io.File
import java.io.PrintWriter

/**
 * Reads only SNPs from drug resistance list, and produces markers
 *
 *
 */

object DrugResistanceMarkers {
  val usage = "scala DrugResistanceMarkers.scala [Drug resistance list] [reference fasta] [marker output]"
  def main(args: Array[String]) {

    if (args.size != 3) println("Wrong argument")
    else {

      val ref = Source.fromFile(args(1)).getLines.filterNot(_.startsWith(">")).mkString

      class SNP(val drug: String, val r: String, val p: Int, val a: String) extends Ordered[SNP] {
        override def toString(): String = ">" + r + p + a + "_resistance_" + drug
        def complementToString(drugs: String): String = ">" + r + p + r + "_susceptibility_" + drugs
        def compare(that: SNP): Int = this.p compare that.p
        val marker = ref.substring(p - 11, p - 1) + a + ref.substring(p, p + 10)
        //val cMarker = ref.substring(p - 11, p + 10)
      }

      object SNP {
        def unapply(s: String): Option[(String, String, Int, String)] = {
          val sArr = s.mkString.split("\t")
          val chrPos = sArr(2)
          if (!chrPos.contains("/") && !chrPos.equals("-")) {
            val drug = sArr(0)
            val nChange = sArr(4)
            val ncArr = nChange.split("/")
            val r = ncArr(0)
            val a = ncArr(1)
            val nucleotides = Array[String]("A", "C", "T", "G")
            if (nucleotides.contains(r) && nucleotides.contains(a))
              Some((drug, r, chrPos.toInt, a))
            else None
          } else None
        }
      }

      class Mutation(s: String) {
        def isSNP: Boolean = s match {
          case SNP(d, r, p, a) => true
          case _ => false
        }
      }

      implicit def seqtoBool(s: String) = new Mutation(s)

      val drList = Source.fromFile(new File(args(0))).getLines.filterNot(_.startsWith("#")).filter(_.isSNP).map { line =>
        line match {
          case SNP(d, r, p, a) => new SNP(d, r, p, a)
        }
      }.toList

      
      /**
       * SNPs grouped by position and mutation.
       */

      val snpsPerPosition = drList.groupBy(_.p).map(snp => snp match {
        case (p, snps) => (p, snps.groupBy(_.a))
      }).toList.sortBy(_._1)

      println(snpsPerPosition.flatMap(_._2.flatMap(_._2)).size + " SNPs")

      
      /**
       * Unique SNPs per position. Duplicate mutations indicating resistance to multiple drugs are merged.
       */

      val uniqueSnps = snpsPerPosition.map(pos => pos match {
        case (p, snpMap) => (p, snpMap.map(mutation => mutation match {
          case (alt, snps) => snps.size match {
            case x if x > 1 =>
              val drugs = snps.map(_.drug).mkString("_"); (alt, List(new SNP(drugs, snps.head.r, p, alt)))
            case 1 => (alt, snps)
          }
        }).flatMap(_._2).toList)
      })

      //uniqueSnps.foreach(println)
      println(uniqueSnps.flatMap(_._2).size + " unique SNPs")

      
      /**
       * Convert unique SNPs to markers, and add susceptible marker per genome position.
       */
      
      val markers = uniqueSnps.flatMap(_ match {
        case (pos, snps) =>
          val mList = snps.map(snp => (snp.toString, snp.marker))
          val drugs = snps.flatMap(d => if (d.drug.contains("_")) d.drug.split("_").toList else List(d.drug)).distinct.mkString("_")
          mList :+ (snps.head.complementToString(drugs), ref.substring(pos - 11, pos + 10))
      })

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

      /**
       * Print unique markers to output file.
       */

      
      val pw = new PrintWriter(new File(args(2)))
      pw.println("# Drug resistance markers")
      markers.foreach(_ match {
        case (id, marker) => pw.println(id + "\n" + marker)
      })
      pw.close

    }

  }
}