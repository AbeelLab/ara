package ara

import java.io.File
import scala.io.Source
import java.io.PrintWriter
import ara.Mutation._

object MutationRate {
  
  case class Config(
    val pathFile: File = null,
    val outFile : File = null
  )
  
  
  def main(args: Array[String]) {
    
    val parser = new scopt.OptionParser[Config]("java -jar ara.jar mutation-rate"){
      opt[File]('i', "input") required() action { (x, c) => c.copy(pathFile = x) } text ("Input file with list of directory paths.")
      opt[File]('o', "output") required() action { (x, c) => c.copy(outFile = x) } text ("Output File to save results.")
      
    }
    
    parser.parse(args, Config()) map { config => 
            
      val vcfs = Source.fromFile(config.pathFile).getLines.map(new File(_)).toList
      val pw = new PrintWriter(config.outFile)
      val drMarkers = scala.io.Source.fromInputStream(MutationRate.getClass().getResourceAsStream("/drMarkers.txt")).getLines().filterNot(f => f.startsWith("#") || f.trim.size == 0).grouped(2).toList
      
      
      val mutations: List[Mutation] = vcfs.map{file =>
        val mutations = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).filter(_.isValid).map(_ match {
          case SNP(r, c, a) => new SNP(r, c, a) 
          case Insertion(r, c, a) => new Insertion(r, c, a) 
          case Deletion(r, c, a) => new Deletion(r, c, a) 
          case MNP(r, c, a) => new MNP(r, c, a)
        })
        mutations
      }.toList.flatten
      
      println(mutations.size + " valid mutations")
      
      val snps = mutations.flatMap{
        case m: SNP => Some(m)
        case _ => None
      }
      println(snps.size + " snps.")
      
      val insertions = mutations.flatMap{
        case m: Insertion => Some(m)
        case _ => None
      }
      println(insertions.size + " insertions.")
      
      val deletions = mutations.flatMap{
        case m: Deletion => Some(m)
        case _ => None
      }
      println(deletions.size + " deletions.")
      
      val mnps = mutations.flatMap{
        case m: MNP => Some(m)
        case _ => None
      }
      println(deletions.size + " deletions.")
      
      val windows = drMarkers.map(s => s(0).split("_")(0).drop(2).dropRight(1).toInt).distinct
      println("drug positions: " + windows.size)
      
      
      val windowCounts = windows.map(i => (i, 
          snps.count(m => Range(i - 10, i + 11).contains(m.coordination)), 
          deletions.count(m => Range(i - 10, i + 11).contains(m.coordination)), 
          insertions.count(m => Range(i - 10, i + 11).contains(m.coordination)),
          mnps.count(m => Range(i - 10, i + 11).contains(m.coordination))
      ))
      
      pw.println("# DR SNP position\tSNPs\tInsertions\tDeletions\tMNPs")
      windowCounts.foreach(m => pw.println(m._1 + "\t" + m._2 + "\t" + m._3 + "\t" + m._4 + "\t" + m._5))
      
      pw.close
    }
    
    
  }
}