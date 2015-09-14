package macaw2

import scala.io.Source
import java.io.File

/**
 * Check whether there are snp-sequences in a phy-file with more than 5% of N's in the total sequence. 
 * 
 */

object QCsnpSequences {
  def main(args: Array[String]) {

    def countN(seq: String): Int = {
      seq.filter(p => p == 'N').size
    }

    if (args.length != 1)
      println("usage: java -jar ara.jar qc-snp-sequences [phy-file]")
    else {
      val phyFile = Source.fromFile(args(0))
      val lineIterator = phyFile.getLines
      val totalSeqSize = lineIterator.take(1).mkString.split(" ")(1).toInt
      val seqIterator = lineIterator.drop(1)
      val countList = lineIterator.map(line => (line.substring(0, 9), countN(line.substring(10)))).toList
      
      val excludeList = countList.filter(s => (s._2.toFloat / totalSeqSize) >= 0.05)
      println(excludeList.size + " Files to exclude.")  
      excludeList.foreach(s => println(s._1 + " with " + s._2 + " Ns.") )
    }

  }
}