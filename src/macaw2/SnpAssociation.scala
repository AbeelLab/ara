//package macaw2

import scala.io.Source
import java.io.File
import java.io.PrintWriter

/**
 * @author Arlin
 */
object SnpAssociation {

  val usage = "scala SnpAssociation.scala [clusterfile] [pathfile] [reference fastafile]"

  def main(args: Array[String]) {

    /** Elapsed time function */
    def time[R](block: => R): R = {
      val t0 = System.currentTimeMillis()
      val result = block // call-by-name
      val t1 = System.currentTimeMillis()
      println("Elapsed time: " + (t1 - t0) + "ms")
      result
    }

    if (args.size != 3) println(usage) else time {

      object SNP {
        def unapply(s: String): Option[(String, Int, String)] = {
          val arr = s.mkString.split("\t")
          val r = arr(3)
          val a = arr(4)
          if (arr(6) == "PASS" && r.length() == 1 && a.length() == 1)
            Some((r, arr(1).toInt, a))
          else None
        }
      }

      def isSNP(line: String): Boolean = line match {
        case SNP(r, p, a) => true
        case _ => false
      }

      def invalidSite(line: String): Boolean = {
        val arr = line.split("\t")
        if (arr(4) == "." && arr(6) == "PASS") true
        else false
      }
      
      
      /**
       * Map of clusters. Key = cluster name, value = List of sample names
       */
      val it = Source.fromFile(args(0)).getLines
      val ls = it.map { line => val l = line.split("\t"); (l(0), l(1)) }.toList
      val clusters = ls.groupBy(f => f._2).mapValues(v => v.map(s => s._1).toList)

      val fileList = Source.fromFile(new File(args(1))).getLines.map(new File(_)).toList

      /**
       * Map with all samples and their list of SNPs
       */
      val snpLists = fileList.map { file =>
        val name = file.getParentFile.getName
        val snps = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filterNot(invalidSite(_)).filter(isSNP(_)).map(_ match {
          case SNP(r, p, a) => (r, p, a)
        }).toList
        (name, snps)
      }.toMap

      val totalSnps = snpLists.flatMap(s => s._2).toList
      val totalDistinctSnps = totalSnps.distinct.sortBy(snp => snp._2)
      println(snpLists.size + " Samples, " + totalSnps.size + " SNPs in total, of which " + totalDistinctSnps.size + " distinct SNPs.")      

      clusters.map{ c =>
        c match {
          case (cName, cList) => {
            val c = cList.filterNot(_ == "MT_H37RV_BRD_V5").flatMap(sample => snpLists(sample))
            println(cName + ": " + c.size + " SNPs, of which " + c.distinct.size + " distinct SNPs.")
          }
        }
      }
      
      
      /**
       * Sublineage association
       */
      val associatedSnps = clusters.map { c =>
        c match {
          case (cName, cList) => {
            val cSnpCounts = cList.filterNot(_ == "MT_H37RV_BRD_V5").flatMap(sample => snpLists(sample)).groupBy(identity).mapValues(_.size) // All SNPs in this cluster and their counts
            val notcNames = (clusters - cName).flatMap(s => s._2).toList // All other clusters names.
            val notcSnpCounts = notcNames.filterNot(_ == "MT_H37RV_BRD_V5").flatMap(sample => snpLists(sample)).groupBy(identity).mapValues(_.size)
            if (cList.contains("MT_H37RV_BRD_V5")) { //Inverse SNPs indicating the absence of this cluster.
              val notcSNPs95 = notcSnpCounts.filter(s => s match { case (snp, count) => (count > notcNames.size * 0.95) }).map(kv => kv._1).toList // SNPs in more than 95% of samples in all other clusters.
              (cName, notcSNPs95.map(snp => if (cSnpCounts.contains(snp)) snp -> cSnpCounts(snp) else snp -> 0).filter(s => s match {case (snp, count) => (count < 0.05 * cName.size)}).map(kv => kv._1))
            } else { // SNPs indicating the presence of this cluster.
              val cSNPs95 = cSnpCounts.filter(s => s match { case (snp, count) => (count > cList.size * 0.95) }).map(kv => kv._1).toList  // SNPs in more than 95% of samples in cluster.
              (cName, cSNPs95.map(snp => if (notcSnpCounts.contains(snp)) snp -> notcSnpCounts(snp) else snp -> 0).filter(s => s match {case (snp, count) => (count < 0.05 * notcNames.size)}).map(kv => kv._1))
            }
          }
        }
      }

      associatedSnps.toList.sortBy(_._1).foreach{ c =>
        c match {
          case (cName, snpList) => {
            println(cName + ": " + snpList.size + " SNPs specific to this cluster.")
          }
        }        
      }
      
      val associatedSnpsPos = associatedSnps.flatMap(c => c._2).map(_._2).toList.sorted
      //println(associatedSnpsPos)
      println(associatedSnpsPos.size + " cluster specific SNPs.")

      /**
       * Remove SNP positions within 10 bp
       */
      val associatedSnpsPos2 = associatedSnpsPos.filterNot { x => 
        val idx = associatedSnpsPos.indexOf(x)
        if (idx == 0) (associatedSnpsPos(idx + 1) - associatedSnpsPos(idx) < 11)
        else if (idx == associatedSnpsPos.size - 1) (associatedSnpsPos(idx) - associatedSnpsPos(idx - 1) < 11)
        else (associatedSnpsPos(idx) - associatedSnpsPos(idx - 1) < 11) || (associatedSnpsPos(idx + 1) - associatedSnpsPos(idx) < 11)
      }
      
      val associatedSnps2 = associatedSnps.map{c => 
        c match {
          case (cName, snpList) => {
            (cName -> snpList.filter(snp => associatedSnpsPos2.contains(snp._2)))
          }
        }  
      }
      
      associatedSnps2.toList.sortBy(_._1).foreach{ c =>
        c match {
          case (cName, snpList) => {
            println(cName + ": " + snpList.size + " SNPs not within 10 bp.")
          }
        }        
      }
      //println(associatedSnpsPos2)
      println(associatedSnpsPos2.size + " non-overlapping SNPs.")

      /**
       * Generate markers
       */
      
      val ref = Source.fromFile(args(2)).getLines.filterNot(_.startsWith(">")).mkString //Source.fromInputStream(getClass.getResourceAsStream("MT_H37RV_BRD_V5.fasta"))
      val markers = associatedSnps2.map { c =>
        c match {
          case (cName, snpList) => {
            if (clusters(cName).contains("MT_H37RV_BRD_V5")) {
              (cName -> snpList.map(snp => (">" + snp._1 + snp._2 + snp._3 + "_absence_" + cName, ref.substring(snp._2 - 11, snp._2 - 1) + snp._3 + ref.substring(snp._2, snp._2 + 10))))
            } else {
              (cName -> snpList.map(snp => (">" + snp._1 + snp._2 + snp._3 + "_presence_" + cName, ref.substring(snp._2 - 11, snp._2 - 1) + snp._3 + ref.substring(snp._2, snp._2 + 10))))
            }
          }
        }
      }

      /**
       * Remove non-unique markers
       */

      val allMarkers = markers.flatMap(c => c._2).map(_._2).toList
      val mCounts = allMarkers.map(m1 => (m1, allMarkers.count(m2 => m2 == m1))).toMap
      //println("Markers + counts: " + mCounts)
      
      val selection = markers.map{ c =>
        c match {
          case (cName, snpList) => {
            (cName, snpList.filter(snp => mCounts(snp._2) == 1))
          }
        }  
      }

      selection.toList.sortBy(_._1).foreach{c => 
        c match {
          case (cName, snpList) => {
            println(cName + ": " + snpList.size + " unique markers.")
          }
        }
      }
      
      println(selection.flatMap(_._2).size + " unique markers in total.")

      /**
       * Print SNP selection to file
       */
      val pw = new java.io.PrintWriter(new File("SNPinfo.txt"))
      pw.println("# MTBC sublineage markers")
      selection foreach { c =>
        c match {
          case (cName, mList) => {
            mList.foreach(m => pw.println(m._1 + "\n" + m._2))
          }
        }
      } 
      pw.close

    }
  }
}