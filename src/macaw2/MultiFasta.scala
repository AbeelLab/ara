package macaw2

import java.io.File
import java.io.PrintWriter
import scala.io.Source

/**
 * @author Arlin
 */
object MultiFasta {
  val usage = "scala MultiFasta.scala [directory] [output.fa] \nMultiple directories may be given as input and should be separated with spaces."

  def main(args: Array[String]) {
    /**
     * Function to write in new textfile.
     */
    def printToFile(f: File)(op: PrintWriter => Unit) {
      val writer = new PrintWriter(f)
      try { op(writer) } finally { writer.close() }
    }

    /**
     * List all Fasta files in the given directory.
     */
    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles.toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.equals("pilon.fasta")) => List(f)
      case _ => Nil
    }

    if (args.length < 1) {
      println(usage)
    } else {
      val list = (0 to args.length - 2).toList.flatMap(arg => listFiles(new File(args(arg))))
      printToFile(new File(args(args.length - 1))) { p =>
        list.foreach { f =>
          p.println(">" + f.getParentFile.getName)
          val it = Source.fromFile(f).getLines
          it.filterNot(line => line.startsWith(">")).foreach(line => p.println(line))
        }
      }
      println("Total files of " + list.length + " fasta-files converted into a multi fasta-file named " + args(args.length - 1) + ".")
    }

  }
}